package nl.esciencecenter.asterisk;

import java.util.concurrent.TimeUnit;

import javax.swing.JFrame;

import nl.esciencecenter.asterisk.data.GlueTimedPlayer;
import nl.esciencecenter.asterisk.input.AsteriskInputHandler;
import nl.esciencecenter.asterisk.interfaces.TimedPlayer;
import nl.esciencecenter.esight.ESightNewtWindow;
import nl.esciencecenter.esight.math.VecF3;
import nl.esciencecenter.esight.math.VecF4;

public class Asterisk {
    private final static AsteriskSettings settings = AsteriskSettings.getInstance();

    private static AsteriskInterfaceWindow amusePanel;
    private static AsteriskGLEventListener amuseWindow;
    private static AsteriskInputHandler amuseInputHandler;

    public Asterisk() {
        // Create the Swing interface elements
        amusePanel = new AsteriskInterfaceWindow();

        // Create the GLEventListener
        amuseWindow = new AsteriskGLEventListener(AsteriskInputHandler.getInstance());

        amuseInputHandler = amuseWindow.getInputHandler();

        new ESightNewtWindow(true, amuseInputHandler, amuseWindow, settings.getDefaultScreenWidth(),
                settings.getDefaultScreenHeight(), "Asterisk - Amuse Visualization Tool");

        // Create the frame
        final JFrame frame = new JFrame("- * -");
        frame.addWindowListener(new java.awt.event.WindowAdapter() {
            @Override
            public void windowClosing(java.awt.event.WindowEvent arg0) {
                System.exit(0);
            }
        });

        frame.setSize(settings.getInterfaceWidth(), settings.getInterfaceHeight());

        frame.setResizable(false);

        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                try {
                    frame.getContentPane().add(amusePanel);
                } catch (final Exception e) {
                    e.printStackTrace(System.err);
                    System.exit(1);
                }
            }
        });

        frame.setVisible(true);
    }

    private static VecF4 randPoint() {
        // Random points on a shpere from http://www.cs.cmu.edu/~mws/rpos.html

        double z = (Math.random() * 2.0) - 1.0;
        double phi = Math.random() * 2.0 * Math.PI;
        double theta = Math.asin(z);
        double ctheta = Math.cos(theta);

        double x = ctheta * Math.cos(phi);
        double y = ctheta * Math.sin(phi);

        VecF4 spherePoint = new VecF4((float) x, (float) y, (float) z, 1f);

        return spherePoint.mul((float) Math.pow(Math.random(), 16));
    }

    private static float randBound(float lower, float delta) {
        return ((float) Math.random() * delta + lower);
    }

    private static float shiftBound(float minBound, float maxBound) {
        float rand = (float) Math.random();
        float x = minBound + (rand * maxBound);

        return (float) Math.pow(x, 3.0);

        // System.err.println(orig + " : " + calc);
    }

    public static void main(String[] args) {
        Asterisk lib = new Asterisk();

        int SCENES = 10;

        int SPHERES = 0;
        int STARS = 100;
        int SPH = 1000;
        int POINTS = 500000;

        float COORD_DELTA_MAX = .3f;
        float COLOR_DELTA_MAX = .2f;
        float RADIUS_DELTA_MAX = .1f;

        // Sphere[] spheres1 = new Sphere[SPHERES];
        // for (int i = 0; i < spheres1.length; i++) {
        // float[] coordinates = new float[] { randBound(-1f, 2f),
        // randBound(-1f, 2f), randBound(-1f, 2f) };
        // float[] color = new float[] { (float) Math.random(),
        // (float) Math.random(), (float) Math.random(), 1f };
        // spheres1[i] = new Sphere(i, coordinates,
        // (float) Math.random() * 0.1f, color);
        // }
        //
        // Star[] stars1 = new Star[STARS];
        // for (int i = 0; i < stars1.length; i++) {
        // float[] coordinates = new float[] { randBound(-1f, 2f),
        // randBound(-1f, 2f), randBound(-1f, 2f) };
        // float[] color = new float[] { (float) Math.random(),
        // (float) Math.random(), (float) Math.random(), 1f };
        // stars1[i] = new Star(i, coordinates, (float) Math.random() * 0.1f,
        // color);
        // }
        //
        // SPHGas[] sphGas1 = new SPHGas[SPH];
        // for (int i = 0; i < sphGas1.length; i++) {
        // float[] coordinates = new float[] { randBound(-1f, 2f),
        // randBound(-1f, 2f), randBound(-1f, 2f) };
        // float[] color = new float[] { (float) Math.random(),
        // (float) Math.random(), (float) Math.random(),
        // (float) Math.random() };
        // sphGas1[i] = new SPHGas(i, coordinates, color);
        // }

        // Perlin3D p3d = new Perlin3D(GL3.GL_TEXTURE0, 100, 100, 100);
        // p3d.init(gl);

        PointGas[] pGas1 = new PointGas[POINTS];
        int tenPercent = (int) (0.1f * pGas1.length);
        long startTime = System.currentTimeMillis();

        for (int i = 0; i < pGas1.length; i++) {
            // float[] coordinates = new float[] { randBound(-1f, 2f),
            // randBound(-1f, 2f), randBound(-1f, 2f) };
            VecF4 rand = randPoint();
            float[] coordinates = new float[] { rand.get(0), rand.get(1), rand.get(2) };
            float[] color = new float[] {
                    // 1f, 1f, 1f, 1f
                    (float) Math.random(), (float) Math.random(), (float) Math.random(), (float) Math.random() };
            pGas1[i] = new PointGas(i, coordinates, color);

            if (i % tenPercent == 0) {
                long stopTime = System.currentTimeMillis();
                long timePassed = stopTime - startTime;

                System.out.println("Point cloud generation progress: "
                        + (int) (((float) i / (float) pGas1.length) * 100) + "% ... Time passed: "
                        + TimeUnit.MILLISECONDS.toSeconds(timePassed) + " seconds.");

            }
        }
        // Snapshot scene1 = new Snapshot("willekeurig", spheres1, stars1,
        // sphGas1, pGas1);

        Snapshot scene1 = new Snapshot("willekeurig", null, null, null, pGas1);

        lib.addScene(scene1);

        // for (int j = 0; j < SCENES; j++) {
        // for (int i = 0; i < spheres1.length; i++) {
        // float[] coordinates = new float[] {
        // shiftBound(spheres1[i].getCoordinates()[0],
        // COORD_DELTA_MAX),
        // shiftBound(spheres1[i].getCoordinates()[1],
        // COORD_DELTA_MAX),
        // shiftBound(spheres1[i].getCoordinates()[2],
        // COORD_DELTA_MAX) };
        //
        // float[] color = new float[] {
        // shiftBound(spheres1[i].getColor()[0], COLOR_DELTA_MAX),
        // shiftBound(spheres1[i].getColor()[1], COLOR_DELTA_MAX),
        // shiftBound(spheres1[i].getColor()[2], COLOR_DELTA_MAX),
        // 1f };
        //
        // spheres1[i] = new Sphere(i, coordinates, shiftBound(
        // spheres1[i].getRadius(), RADIUS_DELTA_MAX), color);
        // }
        //
        // for (int i = 0; i < stars1.length; i++) {
        // float[] coordinates = new float[] {
        // shiftBound(stars1[i].getCoordinates()[0],
        // COORD_DELTA_MAX),
        // shiftBound(stars1[i].getCoordinates()[1],
        // COORD_DELTA_MAX),
        // shiftBound(stars1[i].getCoordinates()[2],
        // COORD_DELTA_MAX) };
        //
        // float[] color = new float[] {
        // shiftBound(stars1[i].getColor()[0], COLOR_DELTA_MAX),
        // shiftBound(stars1[i].getColor()[1], COLOR_DELTA_MAX),
        // shiftBound(stars1[i].getColor()[2], COLOR_DELTA_MAX),
        // 1f };
        //
        // stars1[i] = new Star(i, coordinates, shiftBound(
        // stars1[i].getRadius(), RADIUS_DELTA_MAX), color);
        // }
        //
        // for (int i = 0; i < sphGas1.length; i++) {
        // float[] coordinates = new float[] {
        // shiftBound(sphGas1[i].getCoordinates()[0],
        // COORD_DELTA_MAX),
        // shiftBound(sphGas1[i].getCoordinates()[1],
        // COORD_DELTA_MAX),
        // shiftBound(sphGas1[i].getCoordinates()[2],
        // COORD_DELTA_MAX) };
        //
        // float[] color = new float[] {
        // shiftBound(sphGas1[i].getColor()[0], COLOR_DELTA_MAX),
        // shiftBound(sphGas1[i].getColor()[1], COLOR_DELTA_MAX),
        // shiftBound(sphGas1[i].getColor()[2], COLOR_DELTA_MAX),
        // 1f };
        //
        // sphGas1[i] = new SPHGas(i, coordinates, color);
        // }
        //
        // for (int i = 0; i < pGas1.length; i++) {
        // float[] coordinates = new float[] {
        // shiftBound(pGas1[i].getCoordinates()[0],
        // COORD_DELTA_MAX),
        // shiftBound(pGas1[i].getCoordinates()[1],
        // COORD_DELTA_MAX),
        // shiftBound(pGas1[i].getCoordinates()[2],
        // COORD_DELTA_MAX) };
        //
        // float[] color = new float[] {
        // shiftBound(pGas1[i].getColor()[0], COLOR_DELTA_MAX),
        // shiftBound(pGas1[i].getColor()[1], COLOR_DELTA_MAX),
        // shiftBound(pGas1[i].getColor()[2], COLOR_DELTA_MAX), 1f };
        //
        // pGas1[i] = new PointGas(i, coordinates, color);
        // }
        //
        // Snapshot scene = new Snapshot("random " + j, spheres1, stars1,
        // null, pGas1);
        // lib.addScene(scene);
        // }

        // try {
        // Thread.sleep(3000);
        // makePNGScreenshot("before.png");
        //
        // Thread.sleep(10000);
        // setView(3, 30f, 20f, 0f, -5f);
        // Thread.sleep(3000);
        //
        // makePNGScreenshot("after.png");
        // } catch (InterruptedException e) {
        // // TODO Auto-generated catch block
        // e.printStackTrace();
        // }
    }

    public int addScene(Snapshot scene) {
        TimedPlayer timer = AsteriskInterfaceWindow.getTimer();
        int sceneNumber = ((GlueTimedPlayer) timer).addScene(scene);

        if (!timer.isInitialized()) {
            ((GlueTimedPlayer) timer).init();

            new Thread(timer).start();
        }
        return sceneNumber;
    }

    public void setView(int sceneNumber, float xRotation, float yRotation, float zRotation, float cameraDistance) {
        TimedPlayer timer = AsteriskInterfaceWindow.getTimer();
        ((GlueTimedPlayer) timer).setFrame(sceneNumber, true);

        amuseInputHandler.setRotation(new VecF3(xRotation, yRotation, zRotation));
        amuseInputHandler.setViewDist(cameraDistance);
    }

    public int getCurrentSceneNumberOnDisplay() {
        TimedPlayer timer = AsteriskInterfaceWindow.getTimer();
        return ((GlueTimedPlayer) timer).getFrameNumber();
    }

    public float getXRotation() {
        return amuseInputHandler.getRotation().get(0);
    }

    public float getYRotation() {
        return amuseInputHandler.getRotation().get(1);
    }

    public float getZRotation() {
        return amuseInputHandler.getRotation().get(2);
    }

    public float getCameraDistance() {
        return amuseInputHandler.getViewDist();
    }

    public void makePNGScreenshot(String fileName) {
        // TimedPlayer timer = AsteriskInterfaceWindow.getTimer();
        // ((GlueTimedPlayer) timer).setScreenshotFileName(fileName);
        // ((GlueTimedPlayer) timer).setScreenshotNeeded(true);

        GlueTimedPlayer timer = AsteriskInterfaceWindow.getTimer();
        timer.makeScreenShot(fileName);
    }

    public int getSceneNumber() {
        TimedPlayer timer = AsteriskInterfaceWindow.getTimer();
        return ((GlueTimedPlayer) timer).getFrameNumber();
    }

    public VecF3 getRotation() {
        return amuseInputHandler.getRotation();
    }

    public VecF3 getTranslation() {
        return amuseInputHandler.getTranslation();
    }

    public void setSceneNumber(int sceneNumber) {
        TimedPlayer timer = AsteriskInterfaceWindow.getTimer();
        ((GlueTimedPlayer) timer).setFrame(sceneNumber, true);
    }

    public void setRotation(VecF3 rotation) {
        amuseInputHandler.setRotation(rotation);
    }

    public void setTranslation(VecF3 translation) {
        amuseInputHandler.setTranslation(translation);
    }

    public void setCameraDistance(float cameraDistance) {
        amuseInputHandler.setViewDist(cameraDistance);
    }
}
