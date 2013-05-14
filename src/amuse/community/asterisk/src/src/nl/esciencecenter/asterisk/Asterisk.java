package nl.esciencecenter.asterisk;

import javax.swing.JFrame;

import nl.esciencecenter.asterisk.AsteriskInterfaceWindow;
import nl.esciencecenter.asterisk.data.GlueTimedPlayer;
import nl.esciencecenter.asterisk.input.AsteriskInputHandler;
import nl.esciencecenter.asterisk.interfaces.TimedPlayer;
import nl.esciencecenter.esight.ESightNewtWindow;
import nl.esciencecenter.esight.math.VecF3;

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

    private static float randBound(float lower, float delta) {
        return ((float) Math.random() * delta + lower);
    }

    private static float shiftBound(float orig, float maxDelta) {
        float calc = (float) (orig - (.5f * maxDelta) + (Math.random() * maxDelta));

        // System.err.println(orig + " : " + calc);

        if (calc < -1f) {
            return -1f;
        } else if (calc > 1f) {
            return 1f;
        } else {
            return calc;
        }

    }

    public static void main(String[] args) {
        Asterisk lib = new Asterisk();

        int SCENES = 10;

        int SPHERES = 0;
        int STARS = 100;
        int SPH = 10000;
        int POINTS = 10000;

        float COORD_DELTA_MAX = .3f;
        float COLOR_DELTA_MAX = .2f;
        float RADIUS_DELTA_MAX = .1f;

        Sphere[] spheres1 = new Sphere[SPHERES];
        for (int i = 0; i < spheres1.length; i++) {
            float[] coordinates = new float[] { randBound(-1f, 2f), randBound(-1f, 2f), randBound(-1f, 2f) };
            float[] color = new float[] { (float) Math.random(), (float) Math.random(), (float) Math.random(), 1f };
            spheres1[i] = new Sphere(i, coordinates, (float) Math.random() * 0.1f, color);
        }

        Star[] stars1 = new Star[STARS];
        for (int i = 0; i < stars1.length; i++) {
            float[] coordinates = new float[] { randBound(-1f, 2f), randBound(-1f, 2f), randBound(-1f, 2f) };
            float[] color = new float[] { (float) Math.random(), (float) Math.random(), (float) Math.random(), 1f };
            stars1[i] = new Star(i, coordinates, (float) Math.random() * 0.1f, color);
        }

        SPHGas[] sphGas1 = new SPHGas[SPH];
        for (int i = 0; i < sphGas1.length; i++) {
            float[] coordinates = new float[] { randBound(-1f, 2f), randBound(-1f, 2f), randBound(-1f, 2f) };
            float[] color =
                    new float[] { (float) Math.random(), (float) Math.random(), (float) Math.random(), (float) Math.random() };
            sphGas1[i] = new SPHGas(i, coordinates, color);
        }
        PointGas[] pGas1 = new PointGas[POINTS];
        for (int i = 0; i < pGas1.length; i++) {
            float[] coordinates = new float[] { randBound(-1f, 2f), randBound(-1f, 2f), randBound(-1f, 2f) };
            float[] color = new float[] {
                    // 1f, 1f, 1f, 1f
                    (float) Math.random(), (float) Math.random(), (float) Math.random(), (float) Math.random() };
            pGas1[i] = new PointGas(i, coordinates, color);
        }
        Snapshot scene1 = new Snapshot("willekeurig", spheres1, stars1, sphGas1, pGas1);
        lib.addScene(scene1);

        for (int j = 0; j < SCENES; j++) {
            for (int i = 0; i < spheres1.length; i++) {
                float[] coordinates =
                        new float[] { shiftBound(spheres1[i].getCoordinates()[0], COORD_DELTA_MAX),
                                shiftBound(spheres1[i].getCoordinates()[1], COORD_DELTA_MAX),
                                shiftBound(spheres1[i].getCoordinates()[2], COORD_DELTA_MAX) };

                float[] color =
                        new float[] { shiftBound(spheres1[i].getColor()[0], COLOR_DELTA_MAX),
                                shiftBound(spheres1[i].getColor()[1], COLOR_DELTA_MAX),
                                shiftBound(spheres1[i].getColor()[2], COLOR_DELTA_MAX), 1f };

                spheres1[i] = new Sphere(i, coordinates, shiftBound(spheres1[i].getRadius(), RADIUS_DELTA_MAX), color);
            }

            for (int i = 0; i < stars1.length; i++) {
                float[] coordinates =
                        new float[] { shiftBound(stars1[i].getCoordinates()[0], COORD_DELTA_MAX),
                                shiftBound(stars1[i].getCoordinates()[1], COORD_DELTA_MAX),
                                shiftBound(stars1[i].getCoordinates()[2], COORD_DELTA_MAX) };

                float[] color =
                        new float[] { shiftBound(stars1[i].getColor()[0], COLOR_DELTA_MAX),
                                shiftBound(stars1[i].getColor()[1], COLOR_DELTA_MAX),
                                shiftBound(stars1[i].getColor()[2], COLOR_DELTA_MAX), 1f };

                stars1[i] = new Star(i, coordinates, shiftBound(stars1[i].getRadius(), RADIUS_DELTA_MAX), color);
            }

            for (int i = 0; i < sphGas1.length; i++) {
                float[] coordinates =
                        new float[] { shiftBound(sphGas1[i].getCoordinates()[0], COORD_DELTA_MAX),
                                shiftBound(sphGas1[i].getCoordinates()[1], COORD_DELTA_MAX),
                                shiftBound(sphGas1[i].getCoordinates()[2], COORD_DELTA_MAX) };

                float[] color =
                        new float[] { shiftBound(sphGas1[i].getColor()[0], COLOR_DELTA_MAX),
                                shiftBound(sphGas1[i].getColor()[1], COLOR_DELTA_MAX),
                                shiftBound(sphGas1[i].getColor()[2], COLOR_DELTA_MAX), 1f };

                sphGas1[i] = new SPHGas(i, coordinates, color);
            }

            for (int i = 0; i < pGas1.length; i++) {
                float[] coordinates =
                        new float[] { shiftBound(pGas1[i].getCoordinates()[0], COORD_DELTA_MAX),
                                shiftBound(pGas1[i].getCoordinates()[1], COORD_DELTA_MAX),
                                shiftBound(pGas1[i].getCoordinates()[2], COORD_DELTA_MAX) };

                float[] color =
                        new float[] { shiftBound(pGas1[i].getColor()[0], COLOR_DELTA_MAX),
                                shiftBound(pGas1[i].getColor()[1], COLOR_DELTA_MAX),
                                shiftBound(pGas1[i].getColor()[2], COLOR_DELTA_MAX), 1f };

                pGas1[i] = new PointGas(i, coordinates, color);
            }

            Snapshot scene = new Snapshot("random " + j, spheres1, stars1, null, pGas1);
            lib.addScene(scene);
        }

        try {
            Thread.sleep(3000);
            lib.makePNGScreenshot("before.png");

            Thread.sleep(10000);
            
            lib.setSceneNumber(3);
            lib.setRotation(30f, 20f, 0f);
            lib.setCameraDistance(-5f);
            Thread.sleep(3000);

            lib.makePNGScreenshot("after.png");
        } catch (InterruptedException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    public void addScene(Snapshot scene) {
        TimedPlayer timer = AsteriskInterfaceWindow.getTimer();
        ((GlueTimedPlayer) timer).addScene(scene);

        if (!timer.isInitialized()) {
            ((GlueTimedPlayer) timer).init();

            new Thread(timer).start();
        }

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

    public float getCameraDistance() {
        return amuseInputHandler.getViewDist();
    }

    public void makePNGScreenshot(String fileName) {
        TimedPlayer timer = AsteriskInterfaceWindow.getTimer();
        
        timer.makeScreenShot(fileName);
        
//        ((GlueTimedPlayer) timer).setScreenshotFileName(fileName);
//        ((GlueTimedPlayer) timer).setScreenshotNeeded(true);
    }

    //new

    public void setSceneNumber(int sceneNumber) {
        TimedPlayer timer = AsteriskInterfaceWindow.getTimer();
        ((GlueTimedPlayer) timer).setFrame(sceneNumber, true);
    }

    public void setRotation(float x, float y, float z) {
        amuseInputHandler.setRotation(new VecF3(x, y, z));
    }

    public void setTranslation(float x, float y, float z) {
        amuseInputHandler.setTranslation(new VecF3(x, y, z));
    }

    public void setCameraDistance(float cameraDistance) {
        amuseInputHandler.setViewDist(cameraDistance);
    }

}
