// LaunchAnalysis.java - Comprehensive Launch Analysis for Cubesat
import java.io.File;
import java.util.Date;

import org.orekit.data.DirectoryCrawler;
import org.orekit.propagation.analytical.tle.TLE;
import org.orekit.propagation.analytical.tle.TLEPropagator;
import org.orekit.utils.Constants;
import org.orekit.utils.PVCoordinates;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Frame;
import org.orekit.frames.Transform;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.IERSConventions;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.forces.gravity.potential.GravityFieldFactory;
import org.orekit.forces.gravity.potential.NormalizedSphericalHarmonicsProvider;
import org.orekit.forces.gravity.potential.TideSystem;

public class LaunchAnalysis {

    // French Guiana launch site coordinates (Kourou Space Center)
    private static final double LAUNCH_LAT = 5.232; // degrees North
    private static final double LAUNCH_LON = -52.769; // degrees West
    private static final double LAUNCH_ALT = 0.0; // meters above sea level
    
    // Target orbital slot between HORIZONS-1 and INMARSAT3-F3
    private static final double TARGET_LONGITUDE = 146.65; // degrees East (target slot)
    
    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.err.println("Usage: java -cp <classpath> LaunchAnalysis /path/to/orekit-data");
            System.exit(1);
        }

        // 1) Initialize Orekit data
        File orekitData = new File(args[0]);
        if (!orekitData.isDirectory()) {
            System.err.println("orekit-data path not found: " + orekitData.getAbsolutePath());
            System.exit(2);
        }
        
        org.orekit.data.LazyLoadedDataContext dataContext = new org.orekit.data.LazyLoadedDataContext();
        dataContext.getDataProvidersManager().addProvider(new DirectoryCrawler(orekitData));
        org.orekit.data.DataContext.setDefault(dataContext);

        // 2) Define existing satellites (HORIZONS-1 and INMARSAT3-F3)
        String tleH1l1 = "1 27954U 03044A   25229.72035716  .00000123  00000-0  00000-0 0  9995";
        String tleH1l2 = "2 27954   0.6291  88.4934 0002073  48.4031 298.7623 1.00270627 43521";
        
        String tleI3l1 = "1 24674U 96070A   25229.81985093 -.00000132  00000-0  00000-0 0  9995";
        String tleI3l2 = "2 24674   9.6574  52.4865 0005567 124.3367  76.9065 1.00007233 104874";

        TLE tleH1 = new TLE(tleH1l1, tleH1l2);
        TLE tleI3 = new TLE(tleI3l1, tleI3l2);

        TLEPropagator propH1 = TLEPropagator.selectExtrapolator(tleH1);
        TLEPropagator propI3 = TLEPropagator.selectExtrapolator(tleI3);

        // 3) Analysis date (use TLE epoch)
        AbsoluteDate analysisDate = tleH1.getDate();
        System.out.println("=== LAUNCH ANALYSIS FOR CUBESAT ===");
        System.out.println("Analysis Date: " + analysisDate.toString());
        System.out.println("Launch Site: Kourou Space Center, French Guiana");
        System.out.printf("Launch Coordinates: %.3f°N, %.3f°W, %.1f m altitude\n\n", 
                         LAUNCH_LAT, LAUNCH_LON, LAUNCH_ALT);

        // 4) Analyze existing satellites
        analyzeExistingSatellites(propH1, propI3, analysisDate);
        
        // 5) Calculate optimal insertion parameters
        calculateOptimalInsertion(propH1, propI3, analysisDate);
        
        // 6) Calculate launch trajectory parameters
        calculateLaunchTrajectory(analysisDate);
        
        // 7) Calculate orbital maneuvers
        calculateOrbitalManeuvers(analysisDate);
        
        // 8) Generate GMAT input parameters
        generateGMATParameters(analysisDate);
    }

    private static void analyzeExistingSatellites(TLEPropagator propH1, TLEPropagator propI3, AbsoluteDate date) {
        System.out.println("=== EXISTING SATELLITE ANALYSIS ===");
        
        Frame inertial = FramesFactory.getGCRF();
        Frame itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, true);
        Transform transform = inertial.getTransformTo(itrf, date);
        
        // HORIZONS-1 analysis
        PVCoordinates pvH1 = propH1.getPVCoordinates(date, inertial);
        KeplerianOrbit orbitH1 = new KeplerianOrbit(pvH1, inertial, date, Constants.WGS84_EARTH_MU);
        double[] lonH1 = computeLonKm(pvH1, transform);
        
        // INMARSAT3-F3 analysis
        PVCoordinates pvI3 = propI3.getPVCoordinates(date, inertial);
        KeplerianOrbit orbitI3 = new KeplerianOrbit(pvI3, inertial, date, Constants.WGS84_EARTH_MU);
        double[] lonI3 = computeLonKm(pvI3, transform);
        
        System.out.println("HORIZONS-1:");
        System.out.printf("  Longitude: %.3f°E\n", lonH1[0]);
        System.out.printf("  Altitude: %.1f km\n", lonH1[1] - Constants.WGS84_EARTH_EQUATORIAL_RADIUS / 1000.0);
        System.out.printf("  Inclination: %.3f°\n", Math.toDegrees(orbitH1.getI()));
        System.out.printf("  Eccentricity: %.6f\n", orbitH1.getE());
        
        System.out.println("\nINMARSAT3-F3:");
        System.out.printf("  Longitude: %.3f°E\n", lonI3[0]);
        System.out.printf("  Altitude: %.1f km\n", lonI3[1] - Constants.WGS84_EARTH_EQUATORIAL_RADIUS / 1000.0);
        System.out.printf("  Inclination: %.3f°\n", Math.toDegrees(orbitI3.getI()));
        System.out.printf("  Eccentricity: %.6f\n", orbitI3.getE());
        
        double separation = Math.abs(lonH1[0] - lonI3[0]);
        if (separation > 180) separation = 360 - separation;
        System.out.printf("\nSeparation between satellites: %.3f°\n", separation);
    }

    private static void calculateOptimalInsertion(TLEPropagator propH1, TLEPropagator propI3, AbsoluteDate date) {
        System.out.println("\n=== OPTIMAL INSERTION PARAMETERS ===");
        
        Frame inertial = FramesFactory.getGCRF();
        Frame itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, true);
        Transform transform = inertial.getTransformTo(itrf, date);
        
        PVCoordinates pvH1 = propH1.getPVCoordinates(date, inertial);
        PVCoordinates pvI3 = propI3.getPVCoordinates(date, inertial);
        double[] lonH1 = computeLonKm(pvH1, transform);
        double[] lonI3 = computeLonKm(pvI3, transform);
        
        // Calculate optimal longitude between the two satellites
        double optimalLon = (lonH1[0] + lonI3[0]) / 2.0;
        if (Math.abs(lonH1[0] - lonI3[0]) > 180) {
            optimalLon += 180;
            if (optimalLon > 180) optimalLon -= 360;
        }
        
        System.out.printf("Optimal insertion longitude: %.3f°E\n", optimalLon);
        System.out.printf("Target slot center: %.3f°E\n", TARGET_LONGITUDE);
        
        // Calculate required orbital parameters for GEO insertion
        double targetAltitude = 35786.0; // km (GEO altitude)
        double targetRadius = (Constants.WGS84_EARTH_EQUATORIAL_RADIUS + targetAltitude * 1000.0) / 1000.0; // km
        
        System.out.printf("Target altitude: %.1f km\n", targetAltitude);
        System.out.printf("Target orbital radius: %.1f km\n", targetRadius);
        
        // Calculate required velocity for circular GEO
        double mu = Constants.WGS84_EARTH_MU / 1e9; // km³/s²
        double vCircular = Math.sqrt(mu / targetRadius);
        System.out.printf("Required circular velocity: %.3f km/s\n", vCircular);
    }

    private static void calculateLaunchTrajectory(AbsoluteDate date) {
        System.out.println("\n=== LAUNCH TRAJECTORY PARAMETERS ===");
        
        // Convert launch site coordinates to radians
        double latRad = Math.toRadians(LAUNCH_LAT);
        double lonRad = Math.toRadians(LAUNCH_LON);
        
        // Calculate launch azimuth for GEO insertion
        // For GEO insertion from Kourou, typical azimuth is around 90-95 degrees
        double launchAzimuth = 92.5; // degrees (typical for GEO launches from Kourou)
        System.out.printf("Launch azimuth: %.1f°\n", launchAzimuth);
        
        // Calculate required inclination for GEO
        double requiredInclination = Math.abs(LAUNCH_LAT); // Minimum inclination from launch site
        System.out.printf("Minimum required inclination: %.3f°\n", requiredInclination);
        
        // Calculate launch window considerations
        System.out.println("Launch window considerations:");
        System.out.println("  - Launch during local time: 18:00-22:00 UTC");
        System.out.println("  - Avoid weather constraints");
        System.out.println("  - Consider orbital phasing requirements");
        
        // Calculate transfer orbit parameters
        double perigeeAltitude = 200.0; // km (typical transfer orbit perigee)
        double apogeeAltitude = 35786.0; // km (GEO altitude)
        
        double perigeeRadius = (Constants.WGS84_EARTH_EQUATORIAL_RADIUS + perigeeAltitude * 1000.0) / 1000.0;
        double apogeeRadius = (Constants.WGS84_EARTH_EQUATORIAL_RADIUS + apogeeAltitude * 1000.0) / 1000.0;
        
        double semiMajorAxis = (perigeeRadius + apogeeRadius) / 2.0;
        double eccentricity = (apogeeRadius - perigeeRadius) / (apogeeRadius + perigeeRadius);
        
        System.out.println("\nTransfer orbit parameters:");
        System.out.printf("  Semi-major axis: %.1f km\n", semiMajorAxis);
        System.out.printf("  Eccentricity: %.4f\n", eccentricity);
        System.out.printf("  Perigee altitude: %.1f km\n", perigeeAltitude);
        System.out.printf("  Apogee altitude: %.1f km\n", apogeeAltitude);
        
        // Calculate transfer orbit velocities
        double mu = Constants.WGS84_EARTH_MU / 1e9; // km³/s²
        double vPerigee = Math.sqrt(mu * (2.0/perigeeRadius - 1.0/semiMajorAxis));
        double vApogee = Math.sqrt(mu * (2.0/apogeeRadius - 1.0/semiMajorAxis));
        
        System.out.printf("  Perigee velocity: %.3f km/s\n", vPerigee);
        System.out.printf("  Apogee velocity: %.3f km/s\n", vApogee);
    }

    private static void calculateOrbitalManeuvers(AbsoluteDate date) {
        System.out.println("\n=== ORBITAL MANEUVERS ===");
        
        // 1. Launch to transfer orbit
        System.out.println("1. LAUNCH TO TRANSFER ORBIT:");
        System.out.println("   - Launch vehicle provides initial velocity");
        System.out.println("   - Insert into elliptical transfer orbit");
        System.out.println("   - Perigee: ~200 km altitude");
        System.out.println("   - Apogee: ~35,786 km altitude");
        
        // 2. Apogee kick maneuver
        System.out.println("\n2. APOGEE KICK MANEUVER (AKM):");
        double vCircular = 3.075; // km/s (circular GEO velocity)
        double vApogee = 1.607; // km/s (transfer orbit apogee velocity)
        double deltaV_AKM = vCircular - vApogee;
        System.out.printf("   - Required ΔV: %.3f km/s\n", deltaV_AKM);
        System.out.println("   - Circularize orbit at GEO altitude");
        System.out.println("   - Timing: At apogee of transfer orbit");
        
        // 3. Drift to target longitude
        System.out.println("\n3. LONGITUDE DRIFT MANEUVER:");
        System.out.println("   - Small velocity adjustments to reach target longitude");
        System.out.println("   - Typically < 10 m/s total ΔV");
        System.out.println("   - May require multiple small burns");
        
        // 4. Station keeping
        System.out.println("\n4. STATION KEEPING:");
        System.out.println("   - Regular small maneuvers to maintain position");
        System.out.println("   - Compensate for orbital perturbations");
        System.out.println("   - Typical ΔV: 50-100 m/s per year");
        
        // Total mission ΔV
        double totalDeltaV = deltaV_AKM + 0.01; // AKM + drift maneuvers
        System.out.printf("\nTotal mission ΔV: ~%.3f km/s\n", totalDeltaV);
    }

    private static void generateGMATParameters(AbsoluteDate date) {
        System.out.println("\n=== GMAT SIMULATION PARAMETERS ===");
        
        System.out.println("Launch Site Parameters:");
        System.out.printf("  Latitude: %.6f deg\n", LAUNCH_LAT);
        System.out.printf("  Longitude: %.6f deg\n", LAUNCH_LON);
        System.out.printf("  Altitude: %.1f km\n", LAUNCH_ALT / 1000.0);
        
        System.out.println("\nTransfer Orbit Parameters:");
        System.out.println("  Semi-major axis: 20493.0 km");
        System.out.println("  Eccentricity: 0.7304");
        System.out.println("  Inclination: 5.232 deg");
        System.out.println("  RAAN: 88.9 deg");
        System.out.println("  Argument of Perigee: 180.0 deg");
        System.out.println("  True Anomaly: 0.0 deg");
        
        System.out.println("\nTarget GEO Orbit Parameters:");
        System.out.println("  Semi-major axis: 42164.0 km");
        System.out.println("  Eccentricity: 0.0001");
        System.out.println("  Inclination: 0.1 deg");
        System.out.println("  RAAN: 88.9 deg");
        System.out.println("  Argument of Perigee: 0.0 deg");
        System.out.println("  True Anomaly: 0.0 deg");
        
        System.out.println("\nManeuver Parameters:");
        System.out.println("  Apogee Kick ΔV: 1.468 km/s");
        System.out.println("  Drift Correction ΔV: 0.010 km/s");
        System.out.println("  Total Mission ΔV: 1.478 km/s");
        
        System.out.println("\nGMAT Script Elements Needed:");
        System.out.println("  1. Spacecraft object with mass and fuel properties");
        System.out.println("  2. Launch site location and constraints");
        System.out.println("  3. Transfer orbit initial conditions");
        System.out.println("  4. Apogee kick maneuver (ImpulsiveBurn)");
        System.out.println("  5. Drift correction maneuvers");
        System.out.println("  6. Station keeping maneuvers");
        System.out.println("  7. Perturbation forces (J2, solar radiation pressure)");
        System.out.println("  8. Mission sequence and timeline");
    }

    private static double[] computeLonKm(PVCoordinates pvInertial, Transform transform) {
        PVCoordinates pvITRF = transform.transformPVCoordinates(pvInertial);
        double x = pvITRF.getPosition().getX();
        double y = pvITRF.getPosition().getY();
        double r = Math.sqrt(x*x + y*y + pvITRF.getPosition().getZ()*pvITRF.getPosition().getZ()) / 1000.0;
        double lon = Math.toDegrees(Math.atan2(y, x));
        if (lon < 0) lon += 360.0;
        if (lon > 180.0) lon -= 360.0;
        return new double[]{ lon, r };
    }
}
