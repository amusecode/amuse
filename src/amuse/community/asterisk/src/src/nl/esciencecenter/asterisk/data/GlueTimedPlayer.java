package nl.esciencecenter.asterisk.data;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import javax.media.opengl.GL3;
import javax.swing.JFormattedTextField;

import nl.esciencecenter.asterisk.AsteriskInterfaceWindow.KeyFrame;
import nl.esciencecenter.asterisk.AsteriskSettings;
import nl.esciencecenter.asterisk.Snapshot;
import nl.esciencecenter.asterisk.input.AsteriskInputHandler;
import nl.esciencecenter.asterisk.interfaces.SceneStorage;
import nl.esciencecenter.asterisk.interfaces.TimedPlayer;
import nl.esciencecenter.esight.math.VecF3;
import nl.esciencecenter.esight.math.VecF4;
import nl.esciencecenter.esight.math.VectorFMath;
import nl.esciencecenter.esight.swing.CustomJSlider;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class GlueTimedPlayer implements TimedPlayer {
    private final GlueDatasetManager dsManager;
    private final GlueSceneStorage sceneStorage;

    private states currentState = states.UNOPENED;
    private final static Logger logger = LoggerFactory.getLogger(GlueTimedPlayer.class);

    private final AsteriskSettings settings = AsteriskSettings.getInstance();
    private int frameNumber;

    private final boolean running = true;
    private boolean initialized = false;
    private boolean fileLessMode = false;

    private long startTime, stopTime;

    private final CustomJSlider timeBar;
    private final JFormattedTextField frameCounter;

    private final AsteriskInputHandler inputHandler;

    private boolean needsScreenshot = false;
    private String screenshotDirectory = "";
    private String screenshotFilename = "";

    private final long waittime = settings.getWaitTimeMovie();

    private class Orientation {
        private final int frameNumber;
        private final VecF3 rotation;
        private final float viewDist;

        public Orientation(int frameNumber, VecF3 rotation, float viewDist) {
            this.frameNumber = frameNumber;
            this.rotation = rotation;
            this.viewDist = viewDist;
        }

        public int getFrameNumber() {
            return frameNumber;
        }

        public VecF3 getRotation() {
            return rotation;
        }

        public float getViewDist() {
            return viewDist;
        }

        @Override
        public String toString() {
            return "#: " + frameNumber + " " + rotation + " " + viewDist;
        }

        @Override
        public int hashCode() {
            int hash = 1;
            hash = hash * 17 + frameNumber;
            return hash;
        }

        @Override
        public boolean equals(Object other) {
            if (other instanceof KeyFrame && ((KeyFrame) other).hashCode() == this.hashCode()) {
                return true;
            } else {
                return false;
            }
        }
    }

    private ArrayList<Orientation> orientationList;

    public GlueTimedPlayer(CustomJSlider timeBar, JFormattedTextField frameCounter) {
        this.timeBar = timeBar;
        this.frameCounter = frameCounter;
        this.inputHandler = AsteriskInputHandler.getInstance();
        this.dsManager = new GlueDatasetManager(1, 4);
        this.sceneStorage = dsManager.getSceneStorage();
    }

    @Override
    public void init() {
        this.fileLessMode = true;
        initialized = true;
    }

    @Override
    public void run() {
        if (!initialized && !fileLessMode) {
            System.err.println("HDFTimer started while not initialized.");
            System.exit(1);
        }

        inputHandler.setRotation(new VecF3(settings.getInitialRotationX(), settings.getInitialRotationY(), 0f));
        inputHandler.setViewDist(settings.getInitialZoom());

        int frame = settings.getInitialSimulationFrame();
        updateFrame(frame, true);

        stop();

        while (running) {
            setScreenshotFileName(frameNumber + ".png");

            if ((currentState == states.PLAYING) || (currentState == states.REDRAWING)
                    || (currentState == states.MOVIEMAKING) || (currentState == states.REVIEW)) {
                try {
                    if (!isScreenshotNeeded()) {
                        startTime = System.currentTimeMillis();

                        if (currentState == states.MOVIEMAKING || currentState == states.REVIEW) {
                            for (Orientation o : orientationList) {
                                if (o.getFrameNumber() == frameNumber) {
                                    VecF3 rotation = o.getRotation().clone();
                                    float viewDist = o.getViewDist();
                                    inputHandler.setRotation(rotation);
                                    inputHandler.setViewDist(viewDist);

                                }
                                if (currentState == states.MOVIEMAKING) {
                                    setScreenshotFileName(frameNumber + ".png");
                                    setScreenshotNeeded(true);
                                }
                            }
                        }

                        // Forward frame
                        if (currentState != states.REDRAWING) {
                            int newFrameNumber;
                            try {
                                newFrameNumber = dsManager.getNextFrameNumber(frameNumber);
                                if (sceneStorage.doneWithLastRequest()) {
                                    updateFrame(newFrameNumber, false);
                                }
                            } catch (IOException e) {
                                // currentState = states.WAITINGONFRAME;
                                // stop();
                                logger.debug("nextFrame returned IOException.");
                            }
                        }

                        // Wait for the _rest_ of the timeframe
                        stopTime = System.currentTimeMillis();
                        long spentTime = stopTime - startTime;

                        if (spentTime < waittime) {
                            Thread.sleep(waittime - spentTime);
                        }
                    }
                } catch (final InterruptedException e) {
                    System.err.println("Interrupted while playing.");
                }
            } else if (currentState == states.STOPPED) {
                try {
                    Thread.sleep(100);
                } catch (final InterruptedException e) {
                    System.err.println("Interrupted while stopped.");
                }
            } else if (currentState == states.REDRAWING) {
                currentState = states.STOPPED;
            } else if (currentState == states.WAITINGONFRAME) {
                try {
                    Thread.sleep(10);
                } catch (final InterruptedException e) {
                    System.err.println("Interrupted while waiting.");
                }
            }
        }
    }

    public synchronized int addScene(Snapshot scene) {
        int sceneNumber = dsManager.addScene(scene);

        timeBar.setMaximum(dsManager.getNumFrames() - 1);
        timeBar.invalidate();

        if (currentState == states.WAITINGONFRAME) {
            start();
        }

        return sceneNumber;
    }

    @Override
    public synchronized void setFrame(int value, boolean overrideUpdate) {
        stop();

        try {
            updateFrame(dsManager.getFrameNumberOfIndex(value), overrideUpdate);
        } catch (IndexNotAvailableException e) {
            e.printStackTrace();
        }
    }

    protected synchronized void updateFrame(int newFrameNumber, boolean overrideUpdate) {
        if (dsManager != null) {
            if (newFrameNumber != frameNumber || overrideUpdate) {

                frameNumber = newFrameNumber;
                settings.setCurrentFrameNumber(newFrameNumber);

                this.timeBar.setValue(dsManager.getIndexOfFrameNumber(newFrameNumber));
                this.frameCounter.setValue(dsManager.getIndexOfFrameNumber(newFrameNumber));
            }
        }
    }

    @Override
    public synchronized void oneBack() {
        stop();

        try {
            int newFrameNumber = dsManager.getPreviousFrameNumber(frameNumber);
            updateFrame(newFrameNumber, false);
        } catch (IOException e) {
            logger.debug("one back failed.");
        }
    }

    @Override
    public synchronized void oneForward() {
        stop();

        try {
            int newFrameNumber = dsManager.getNextFrameNumber(frameNumber);
            updateFrame(newFrameNumber, false);
        } catch (IOException e) {
            logger.debug("one forward failed.");
        }
    }

    @Override
    public synchronized void redraw() {
        if (initialized) {
            updateFrame(frameNumber, true);
            currentState = states.REDRAWING;
        }
    }

    @Override
    public synchronized void rewind() {
        stop();
        updateFrame(0, false);
    }

    @Override
    public SceneStorage getSceneStorage() {
        return sceneStorage;
    }

    @Override
    public synchronized void start() {
        currentState = states.PLAYING;
    }

    @Override
    public synchronized void stop() {
        currentState = states.STOPPED;
    }

    @Override
    public boolean isInitialized() {
        return initialized;
    }

    @Override
    public synchronized boolean isPlaying() {
        if ((currentState == states.PLAYING) || (currentState == states.MOVIEMAKING) || (currentState == states.REVIEW)) {
            return true;
        }

        return false;
    }

    @Override
    public synchronized void movieMode() {
        currentState = states.MOVIEMAKING;
    }

    @Override
    public synchronized void reviewMode() {
        currentState = states.REVIEW;
    }

    @Override
    public synchronized void setScreenshotNeeded(boolean value) {
        needsScreenshot = value;
        notifyAll();
    }

    @Override
    public synchronized boolean isScreenshotNeeded() {
        return needsScreenshot;
    }

    public synchronized void setScreenshotFileName(String screenshotFilename) {
        this.screenshotFilename = screenshotDirectory + screenshotFilename;
    }

    @Override
    public synchronized String getScreenshotFileName() {
        return screenshotFilename;
    }

    @Override
    public void close() {
        initialized = false;
        frameNumber = 0;
        timeBar.setValue(0);
        frameCounter.setValue(0);
        timeBar.setMaximum(0);
    }

    @Override
    public int getFrameNumber() {
        return frameNumber;
    }

    @Override
    public void delete(GL3 gl) {
        // TODO Auto-generated method stub
    }

    @Override
    public synchronized void makeScreenShot(String screenshotFilename) {
        this.screenshotFilename = screenshotDirectory + screenshotFilename;
        this.needsScreenshot = true;

        while (this.needsScreenshot) {
            try {
                wait();
            } catch (InterruptedException e) {
                // IGNORE
            }
        }

    }

    @Override
    public void startSequence(ArrayList<KeyFrame> keyFrames, boolean record) {
        int startKeyFrameNumber = 0, finalKeyFrameNumber = 0;
        for (KeyFrame keyFrame : keyFrames) {
            int frameNumber = keyFrame.getFrameNumber();
            if (frameNumber > finalKeyFrameNumber) {
                finalKeyFrameNumber = frameNumber;
            }
            if (frameNumber < startKeyFrameNumber) {
                startKeyFrameNumber = frameNumber;
            }
        }

        orientationList = new ArrayList<Orientation>();

        ArrayList<Integer> intermediateFrameNumbers = new ArrayList<Integer>();

        try {
            int currentFrameNumber = startKeyFrameNumber;
            while (currentFrameNumber <= finalKeyFrameNumber) {
                intermediateFrameNumbers.add(currentFrameNumber);
                currentFrameNumber = dsManager.getNextFrameNumber(currentFrameNumber);
            }
        } catch (IOException e) {
            // We're done.
        }

        // Only do interpolation step if we have a current AND a next keyFrame
        // available.
        if (keyFrames.size() > 1) {
            for (int i = 0; i < keyFrames.size() - 1; i++) {
                KeyFrame currentKeyFrame = keyFrames.get(i);
                KeyFrame nextKeyFrame = keyFrames.get(i + 1);

                int startFrameNumber = currentKeyFrame.getFrameNumber();
                int stopFrameNumber = nextKeyFrame.getFrameNumber();
                int numberOfInterpolationFrames = 0;
                for (int currentFrameNumber : intermediateFrameNumbers) {
                    if (currentFrameNumber >= startFrameNumber && currentFrameNumber < stopFrameNumber) {
                        numberOfInterpolationFrames++;
                    }
                }

                VecF4 startLocation = new VecF4(currentKeyFrame.getRotation(), 1f);
                VecF4 endLocation = new VecF4(nextKeyFrame.getRotation(), 1f);

                VecF3 startControl = new VecF3();
                VecF3 endControl = new VecF3();

                VecF4[] curveSteps = VectorFMath.bezierCurve(numberOfInterpolationFrames, startLocation, startControl,
                        endControl, endLocation);

                // Patch for zoom
                startLocation = new VecF4(currentKeyFrame.getViewDist(), 0f, 0f, 1f);
                endLocation = new VecF4(nextKeyFrame.getViewDist(), 0f, 0f, 1f);

                VecF4[] zoomSteps = VectorFMath.bezierCurve(numberOfInterpolationFrames, startLocation, startControl,
                        endControl, endLocation);

                for (int j = 0; j < numberOfInterpolationFrames; j++) {
                    int currentFrameNumber = intermediateFrameNumbers.get(currentKeyFrame.getFrameNumber() + j);
                    Orientation newOrientation = new Orientation(currentFrameNumber, curveSteps[j].stripAlpha(),
                            zoomSteps[j].get(0));
                    orientationList.add(newOrientation);
                }
            }
        }

        stop();
        setFrame(startKeyFrameNumber, true);

        if (record) {
            movieMode();
        } else {
            reviewMode();
        }
    }

    public void setScreenshotDirectory(String string) {
        File ssDir = new File(screenshotDirectory);
        ssDir.mkdir();
        this.screenshotDirectory = ssDir.getAbsolutePath();
    }
}
