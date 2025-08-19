import java.io.File;
import java.util.Locale;

import org.orekit.utils.PVCoordinates;

import org.orekit.bodies.BodyShape;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.data.DataContext;
import org.orekit.data.DataProvidersManager;
import org.orekit.data.DirectoryCrawler;
import org.orekit.errors.OrekitException;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Predefined;
import org.orekit.models.earth.Geoid;
import org.orekit.propagation.analytical.KeplerianPropagator;
import org.orekit.propagation.analytical.tle.TLE;
import org.orekit.propagation.analytical.tle.TLEPropagator;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;
import org.orekit.utils.IERSConventions;
import org.orekit.orbits.PositionAngleType;

public class GeoMissionPlanner {

    // ---- CONSTANTS (SI) ----
    static final double MU = 3.986004418e14; // Earth's GM [m^3/s^2]
    static final double R_E = Constants.WGS84_EARTH_EQUATORIAL_RADIUS; // [m]
    static final double GEO_ALT = 35786e3; // GEO altitude above mean equator [m]
    static final double R_GEO = R_E + GEO_ALT;

    public static void main(String[] args) {
        Locale.setDefault(Locale.US);

        try {
            // 0) Set up Orekit data context
            if (args.length > 0) {
                String orekitPath = args[0];
                File orekitData = new File(orekitPath);
                if (orekitData.isDirectory()) {
                    org.orekit.data.LazyLoadedDataContext dataContext = new org.orekit.data.LazyLoadedDataContext();
                    dataContext.getDataProvidersManager().addProvider(new DirectoryCrawler(orekitData));
                    org.orekit.data.DataContext.setDefault(dataContext);
                } else {
                    System.err.println("Orekit data path not found: " + orekitData.getAbsolutePath());
                    System.exit(1);
                }
            }

            // 1) Inputs you can change
            //    Launch site: CSG (Kourou), French Guiana (approx)
            double siteLatDeg = 5.239;          // deg
            double siteLonDeg = -52.768;        // deg (West negative)
            double siteAlt = 10.0;              // meters

            //    Target GEO longitude (slot center) in degrees East; West negative
            double targetLonDeg = -146.65;

            //    Launch epoch (UTC) -> choose something realistic/available
            AbsoluteDate launchEpoch = new AbsoluteDate(2025, 9, 15, 12, 0, 0.0, TimeScalesFactory.getUTC());

            //    Initial GTO after separation (typical “launcher-delivered” GTO):
            double perigeeAlt = 250e3;             // 250 km
            double apogeeAlt = GEO_ALT;            // apogee at GEO altitude
            double incGTOdeg = 6.0;                // ~5-7 deg typical from CSG
            //    RAAN/AOP/TA: if you don’t have phasing from launcher, start simple:
            double raanDeg = 0.0;
            double argPerDeg = 0.0;
            double trueAnomDeg = 0.0;

            //    Desired final GEO inclination (stationkeeping target)
            double finalIncDeg = 0.05;

            //    Whether to compute the perigee injection from a 200 km LEO (off-board)
            boolean includePerigeeInjection = false;

            // 2) Neighbors from your TLEs (HORIZONS-1 & INMARSAT-3 F3)
            TLE tleHorizons1 = new TLE(
                "1 27954U 03044A   25229.72035716  .00000123  00000-0  00000-0 0  9995",
                "2 27954   0.6291  88.4934 0002073  48.4031 298.7623  1.00270627 43521"
            );
            TLE tleInmarsat3F3 = new TLE(
                "1 24674U 96070A   25229.81985093 -.00000132  00000-0  00000-0 0  9995",
                "2 24674   9.6574  52.4865 0005567 124.3367  76.9065  1.00007233 104874"
            );

            // Choose a common “assessment” date near launch (or final GEO) to print their longitudes
            AbsoluteDate assessDate = launchEpoch;

            double lonH1 = satSubLonDeg(tleHorizons1, assessDate);
            double lonF3 = satSubLonDeg(tleInmarsat3F3, assessDate);

            // Normalize longitudes to (-180, 180]
            lonH1 = wrap180(lonH1);
            lonF3 = wrap180(lonF3);

            // Gap midpoint (shortest arc)
            double gapMid = wrap180(midLongitude(lonH1, lonF3));

            // 3) Build initial GTO orbit (Keplerian)
            double rp = R_E + perigeeAlt;
            double ra = R_E + apogeeAlt;
            double a = 0.5 * (rp + ra);
            double e = (ra - rp) / (ra + rp);
            double i = Math.toRadians(incGTOdeg);
            double raan = Math.toRadians(raanDeg);
            double argPer = Math.toRadians(argPerDeg);
            double ta = Math.toRadians(trueAnomDeg);

            org.orekit.orbits.KeplerianOrbit gto = new org.orekit.orbits.KeplerianOrbit(
                a, e, i, argPer, raan, ta, PositionAngleType.TRUE,
                FramesFactory.getEME2000(), launchEpoch, MU
            );

            // 4) Time to apogee from current true anomaly (for planning the onboard burn)
            KeplerianPropagator kep = new KeplerianPropagator(gto);
            double E = 2.0 * Math.atan(Math.tan(ta/2.0) * Math.sqrt((1.0 - e)/(1.0 + e)));
            if (E < 0) E += 2*Math.PI;
            double M = E - e * Math.sin(E);
            double n = Math.sqrt(MU / (a*a*a)); // mean motion
            double timeToApogee = (Math.PI - M) / n; // seconds to next apogee (approx)
            AbsoluteDate apogeeDate = launchEpoch.shiftedBy(timeToApogee);

            // 5) Δv at apogee: circularization + plane change to finalIncDeg
            //    Velocity at apogee on GTO:
            double vApogeeGTO = Math.sqrt(MU * (2.0/ra - 1.0/a));
            //    Circular GEO velocity:
            double vGEO = Math.sqrt(MU / R_GEO);
            //    Plane change at apogee from incGTOdeg -> finalIncDeg (assume final RAAN doesn’t matter)
            double dIrad = Math.toRadians(Math.abs(incGTOdeg - finalIncDeg));
            //    Combined burn magnitude:
            double dvApogee = Math.sqrt(vApogeeGTO*vApogeeGTO + vGEO*vGEO - 2.0*vApogeeGTO*vGEO*Math.cos(dIrad));

            double dvPerigee = 0.0;
            if (includePerigeeInjection) {
                // Optional: from a 200 km LEO to GTO perigee burn
                double leoAlt = 200e3;
                double rLEO = R_E + leoAlt;
                double vLEO = Math.sqrt(MU / rLEO);
                double aTrans = a; // same transfer a
                double vPerigeeTrans = Math.sqrt(MU * (2.0/rLEO - 1.0/aTrans));
                dvPerigee = Math.abs(vPerigeeTrans - vLEO);
            }
            double dvTotal = dvPerigee + dvApogee;

            // 6) After apogee burn, you are at near-zero-ecc geo (a ~ R_GEO) with ~finalIncDeg.
            //    Now compute a small drift plan to reach targetLonDeg exactly (if needed).
            //    Simple scheme: set a slightly off-GEO SMA to create a drift rate, then correct.
            double currentSlotDeg = estimateSlotFromGtoApogeePhase(targetLonDeg, lonH1, lonF3, gapMid);
            // but for safety, we just plan to drift from whatever initial longitude at circularization to target
            // Assume we circularize over the apogee longitude of the GTO (rough estimate):
            double lonAtApogee = subLonOfKeplerAtDate(gto, apogeeDate);
            double requiredDeltaLon = wrapSigned(targetLonDeg - lonAtApogee);

            // Drift model: drift rate ≈ (n - nGEO) * 86400 * 180/pi [deg/day], n = sqrt(mu/a^3)
            // For small delta-a around GEO, driftRate ≈ - (3/2) * nGEO * (Δa/a) * 86400 * 180/pi
            double nGEO = Math.sqrt(MU / (R_GEO*R_GEO*R_GEO));
            // Choose a convenient drift rate magnitude (e.g. 0.5 deg/day) and compute needed Δa:
            double desiredDriftDegPerDay = Math.signum(requiredDeltaLon) * 0.5; // adjust sign
            double deltaA = - desiredDriftDegPerDay * (Math.PI/180.0) / (1.5 * nGEO * 86400.0) * R_GEO;

            // Time to drift:
            double daysToTarget = (Math.abs(requiredDeltaLon) / Math.abs(desiredDriftDegPerDay));
            double secondsToTarget = daysToTarget * 86400.0;

            // Small dv to set drift (tangential tweak at GEO):
            // dv ≈ 0.5 * vGEO * |Δa|/a   (first-order)
            double dvSetDrift = 0.5 * vGEO * Math.abs(deltaA) / R_GEO;
            // dv to stop drift at arrival ~ same magnitude:
            double dvStopDrift = dvSetDrift;
            double dvStationkeepingPhase = dvSetDrift + dvStopDrift;

            // 7) Print everything you’ll feed to GMAT
            System.out.println("=== GEO Mission Planner (Orekit) ===");
            System.out.println("Launch site (Kourou): lat=" + siteLatDeg + " deg, lon=" + siteLonDeg + " deg, h=" + siteAlt + " m");
            System.out.println("Launch epoch (UTC):   " + launchEpoch);
            System.out.println();

            System.out.printf("Neighbor longitudes at %s%n", assessDate);
            System.out.printf("  HORIZONS-1 : %7.2f deg%n", lonH1);
            System.out.printf("  INMARSAT-3F3: %7.2f deg%n", lonF3);
            System.out.printf("Gap midpoint : %7.2f deg%n", gapMid);
            System.out.println();

            System.out.println("== Initial GTO (after separation) ==");
            System.out.printf("  SMA a        : %.3f km%n", a/1000.0);
            System.out.printf("  ecc e        : %.6f%n", e);
            System.out.printf("  inc i        : %.3f deg%n", Math.toDegrees(i));
            System.out.printf("  RAAN         : %.3f deg%n", raanDeg);
            System.out.printf("  AOP          : %.3f deg%n", argPerDeg);
            System.out.printf("  TA           : %.3f deg%n", trueAnomDeg);
            System.out.printf("  Perigee alt  : %.1f km%n", perigeeAlt/1000.0);
            System.out.printf("  Apogee  alt  : %.1f km%n", apogeeAlt/1000.0);
            System.out.println();

            System.out.println("== Onboard maneuvers ==");
            if (includePerigeeInjection) {
                System.out.printf("  Δv (perigee injection LEO→GTO): %.2f m/s%n", dvPerigee);
            } else {
                System.out.println("  (Perigee injection provided by launcher; not counted here)");
            }
            System.out.printf("  Time to apogee (from launch): %.1f min%n", timeToApogee/60.0);
            System.out.printf("  Apogee burn epoch (UTC): %s%n", apogeeDate);
            System.out.printf("  Δv (apogee circ + plane change %.2f→%.2f deg): %.2f m/s%n",
                    incGTOdeg, finalIncDeg, dvApogee);
            System.out.printf("  Total Δv (onboard, incl drift set/stop): %.2f m/s%n", (dvApogee + dvStationkeepingPhase));
            System.out.println();

            System.out.println("== GEO drift to slot ==");
            System.out.printf("  GEO circ longitude (rough, at apogee): %7.2f deg%n", lonAtApogee);
            System.out.printf("  Target slot longitude: %7.2f deg%n", targetLonDeg);
            System.out.printf("  Needed Δlon (signed) : %7.2f deg%n", requiredDeltaLon);
            System.out.printf("  Chosen drift rate    : %5.2f deg/day%n", desiredDriftDegPerDay);
            System.out.printf("  Δa for drift         : %+6.0f m (a = Rgeo %+d m)%n", deltaA, 0);
            System.out.printf("  Time to target       : %5.1f days (%.1f h)%n", daysToTarget, daysToTarget*24.0);
            System.out.printf("  dv set/stop drift    : %.2f + %.2f = %.2f m/s%n", dvSetDrift, dvStopDrift, dvStationkeepingPhase);
            System.out.println();

            // 8) Simple separation sanity check: longitudes at final epoch
            AbsoluteDate finalEpoch = apogeeDate.shiftedBy(secondsToTarget);
            double lonH1f = satSubLonDeg(tleHorizons1, finalEpoch);
            double lonF3f = satSubLonDeg(tleInmarsat3F3, finalEpoch);
            double minSepH1 = angularSeparationDeg(wrap180(targetLonDeg), wrap180(lonH1f));
            double minSepF3 = angularSeparationDeg(wrap180(targetLonDeg), wrap180(lonF3f));

            System.out.println("== Separation check at arrival ==");
            System.out.printf("  HORIZONS-1 at %s: %7.2f deg  | separation to target: %.2f deg%n",
                    finalEpoch, wrap180(lonH1f), minSepH1);
            System.out.printf("  INMARSAT-3F3 at %s: %7.2f deg | separation to target: %.2f deg%n",
                    finalEpoch, wrap180(lonF3f), minSepF3);
            System.out.println("  (You’ll still do conjunction screening & ops-level coordination.)");

            // 9) GMAT-ready orbit for the CubeSat (final GEO)
            System.out.println();
            System.out.println("== GMAT final GEO suggestion (CubeSat) ==");
            System.out.printf("  Epoch (UTC): %s%n", finalEpoch);
            System.out.printf("  a (km): %.3f%n", R_GEO/1000.0);
            System.out.printf("  e     : %.5f%n", 0.0001); // near-circular
            System.out.printf("  i (deg): %.3f%n", finalIncDeg);
            System.out.printf("  RAAN (deg): %.3f%n", 0.0); // pick convenient, will not matter in perfect GEO
            System.out.printf("  AOP  (deg): %.3f%n", 0.0);
            System.out.printf("  TA   (deg): %.3f%n", 0.0);
            System.out.printf("  Subsatellite longitude target: %.2f deg%n", targetLonDeg);

        } catch (OrekitException e) {
            System.err.println("Orekit error: " + e.getMessage());
        }
    }

    // ---- Utility: sub-satellite longitude from TLE at date ----
    static double satSubLonDeg(TLE tle, AbsoluteDate date) {
        TLEPropagator prop = TLEPropagator.selectExtrapolator(tle);
        PVCoordinates posITRF = prop.propagate(date).getPVCoordinates(FramesFactory.getITRF(IERSConventions.IERS_2010, true));
        Frame itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, true);
        BodyShape earth = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                                               Constants.WGS84_EARTH_FLATTENING, itrf);
        GeodeticPoint gp = earth.transform(posITRF.getPosition(), itrf, date);
        return Math.toDegrees(gp.getLongitude());
    }

    // ---- Utility: rough sublon of GTO apogee point when circularizing (Kepler-only, no Earth rotation) ----
    static double subLonOfKeplerAtDate(org.orekit.orbits.KeplerianOrbit orb, AbsoluteDate dateAtApogee) {
        // Propagate to apogee using Keplerian model and convert to Earth-fixed longitude
        KeplerianPropagator k = new KeplerianPropagator(orb);
        PVCoordinates posEME = k.propagate(dateAtApogee).getPVCoordinates(FramesFactory.getEME2000());
        // Transform to ITRF (ignoring EOP if data not present)
        Frame itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, true);
        PVCoordinates posITRF = FramesFactory.getEME2000().getTransformTo(itrf, dateAtApogee).transformPVCoordinates(posEME);
        BodyShape earth = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                                               Constants.WGS84_EARTH_FLATTENING, itrf);
        GeodeticPoint gp = earth.transform(posITRF.getPosition(), itrf, dateAtApogee);
        return Math.toDegrees(gp.getLongitude());
    }

    // ---- Longitude helpers ----
    static double wrap180(double lonDeg) {
        double x = lonDeg % 360.0;
        if (x <= -180.0) x += 360.0;
        if (x > 180.0)   x -= 360.0;
        return x;
    }
    static double wrapSigned(double d) {
        double x = d % 360.0;
        if (x > 180) x -= 360;
        if (x < -180) x += 360;
        return x;
    }
    static double midLongitude(double lon1, double lon2) {
        // midpoint along the shortest arc
        double d = wrapSigned(lon2 - lon1);
        return wrap180(lon1 + d/2.0);
    }
    static double angularSeparationDeg(double lonA, double lonB) {
        return Math.abs(wrapSigned(lonA - lonB));
    }

    // Placeholder for slot estimation if you want to auto-snap to gap center
    static double estimateSlotFromGtoApogeePhase(double target, double lon1, double lon2, double gapMid) {
        return gapMid; // keep it simple: aim at gap midpoint
    }
}
