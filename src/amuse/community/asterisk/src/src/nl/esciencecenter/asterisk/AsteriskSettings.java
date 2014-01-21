package nl.esciencecenter.asterisk;

import nl.esciencecenter.asterisk.data.GlueSceneDescription;
import nl.esciencecenter.esight.util.Settings;
import nl.esciencecenter.esight.util.TypedProperties;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class AsteriskSettings extends Settings {
    private final static Logger logger = LoggerFactory.getLogger(AsteriskSettings.class);

    private final double STAR_DEFAULT_LUMINOSITY = 3000.0;

    private static class SingletonHolder {
        public final static AsteriskSettings instance = new AsteriskSettings();
    }

    private final String[] ACCEPTABLE_EXTENSIONS = new String[] { "gas", "bin" };

    // Minimum and maximum values for the brightness sliders
    private float POSTPROCESSING_OVERALL_BRIGHTNESS_MIN = 0f;
    private float POSTPROCESSING_OVERALL_BRIGHTNESS_MAX = 10f;
    private float POSTPROCESSING_AXES_BRIGHTNESS_MIN = 0f;
    private float POSTPROCESSING_AXES_BRIGHTNESS_MAX = 10f;
    private float POSTPROCESSING_POINT_GAS_BRIGHTNESS_MIN = 0f;
    private float POSTPROCESSING_POINT_GAS_BRIGHTNESS_MAX = 10f;
    private float POSTPROCESSING_STAR_HALO_BRIGHTNESS_MIN = 0f;
    private float POSTPROCESSING_STAR_HALO_BRIGHTNESS_MAX = 10f;
    private float POSTPROCESSING_STAR_BRIGHTNESS_MIN = 0f;
    private float POSTPROCESSING_STAR_BRIGHTNESS_MAX = 10f;
    private float POSTPROCESSING_SPHERE_BRIGHTNESS_MIN = 0f;
    private float POSTPROCESSING_SPHERE_BRIGHTNESS_MAX = 10f;
    private float POSTPROCESSING_HUD_BRIGHTNESS_MIN = 0f;
    private float POSTPROCESSING_HUD_BRIGHTNESS_MAX = 10f;

    // Settings for the postprocessing shader
    private float POSTPROCESSING_OVERALL_BRIGHTNESS_DEF = 5f;
    private float POSTPROCESSING_AXES_BRIGHTNESS_DEF = 4f;
    private float POSTPROCESSING_POINT_GAS_BRIGHTNESS_DEF = 10f;

    private float POSTPROCESSING_STAR_HALO_BRIGHTNESS_DEF = 4f;
    private float POSTPROCESSING_STAR_BRIGHTNESS_DEF = 4f;
    private float POSTPROCESSING_SPHERE_BRIGHTNESS_DEF = 4f;
    private float POSTPROCESSING_HUD_BRIGHTNESS_DEF = 4f;

    private final float PARTICLE_SIZE_MULTIPLIER_MIN = 1f;
    private float PARTICLE_SIZE_MULTIPLIER_DEF = 1f;
    private final float PARTICLE_SIZE_MULTIPLIER_MAX = 20f;

    private final int SLIDER_MIN_BLUR_TYPE = 0;
    private final int SLIDER_DEF_BLUR_TYPE = 1;
    private final int SLIDER_MAX_BLUR_TYPE = 9;

    private final int SLIDER_MIN_BLUR_PASSES = 0;
    private final int SLIDER_DEF_BLUR_PASSES = 1;
    private final int SLIDER_MAX_BLUR_PASSES = 8;

    private final int SLIDER_MIN_BLUR_SIZE = 1;
    private final int SLIDER_DEF_BLUR_SIZE = 2;
    private final int SLIDER_MAX_BLUR_SIZE = 10;

    private boolean GAS_COLOR_INVERTED = false;
    private boolean GAS_COLOR_BACKGROUND_INVERTED = false;
    private boolean GAS_COLOR_FROM_STARS = false;

    private boolean STAR_COLORS_EXAGGERATED = true;

    private long WAITTIME_FOR_RETRY = 10000;
    private long WAITTIME_FOR_MOVIE = 200;
    private final float EPSILON = 1.0E-7f;

    private boolean BEZIER_INTERPOLATION = false;
    private int BEZIER_INTERPOLATION_STEPS = 10;

    private int PREPROCESSING_AMOUNT = 5;

    private int pointGasBlurTypeSetting;
    private int pointGasBlurPassSetting;
    private int pointGasBlurSizeSetting;

    private int starHaloBlurTypeSetting;
    private int starHaloBlurPassSetting;
    private int starHaloBlurSizeSetting;

    private float pointGasPointSizeSetting = 1f;

    public static final int STAR_SUBDIVISION = 2;
    public static final int PLANET_SUBDIVISION = 2;
    public static final int SPHERE_SUBDIVISION = 2;

    private float fieldOfView = 45.0f;
    private float zNear = 0.1f;
    private float zFar = 3000.0f;

    private GlueSceneDescription currentDescription;

    public static AsteriskSettings getInstance() {
        return SingletonHolder.instance;
    }

    private final String currentExtension = "bin";
    private final int timeStep = 1;

    private final float SLIDER_MIN_POINT_SIZE = 1f;
    private final float SLIDER_MAX_POINT_SIZE = 10f;

    private boolean pointgasSizeDependantOnCameraDistance = false;

    private boolean doOrbit = false;
    private boolean orbitLinkedToPlayback = false;
    private float orbitSpeed = 0.1f;
    private float xOrbitSpeed = 0f;
    private float yOrbitSpeed = 1f;
    private float zOrbitSpeed = 0f;

    private AsteriskSettings() {
        super();

        try {
            final TypedProperties props = new TypedProperties();
            props.loadFromClassPath("settings.properties");

            POSTPROCESSING_OVERALL_BRIGHTNESS_MIN = props.getFloatProperty("POSTPROCESSING_OVERALL_BRIGHTNESS_MIN");
            POSTPROCESSING_OVERALL_BRIGHTNESS_MAX = props.getFloatProperty("POSTPROCESSING_OVERALL_BRIGHTNESS_MAX");

            // Minimum and maximum values for the brightness sliders
            POSTPROCESSING_OVERALL_BRIGHTNESS_MIN = props.getFloatProperty("POSTPROCESSING_OVERALL_BRIGHTNESS_MIN");
            POSTPROCESSING_OVERALL_BRIGHTNESS_MAX = props.getFloatProperty("POSTPROCESSING_OVERALL_BRIGHTNESS_MAX");
            POSTPROCESSING_AXES_BRIGHTNESS_MIN = props.getFloatProperty("POSTPROCESSING_AXES_BRIGHTNESS_MIN");
            POSTPROCESSING_AXES_BRIGHTNESS_MAX = props.getFloatProperty("POSTPROCESSING_AXES_BRIGHTNESS_MAX");
            POSTPROCESSING_POINT_GAS_BRIGHTNESS_MIN = props.getFloatProperty("POSTPROCESSING_POINT_GAS_BRIGHTNESS_MIN");
            POSTPROCESSING_POINT_GAS_BRIGHTNESS_MAX = props.getFloatProperty("POSTPROCESSING_POINT_GAS_BRIGHTNESS_MAX");
            POSTPROCESSING_STAR_HALO_BRIGHTNESS_MIN = props.getFloatProperty("POSTPROCESSING_STAR_HALO_BRIGHTNESS_MIN");
            POSTPROCESSING_STAR_HALO_BRIGHTNESS_MAX = props.getFloatProperty("POSTPROCESSING_STAR_HALO_BRIGHTNESS_MAX");
            POSTPROCESSING_STAR_BRIGHTNESS_MIN = props.getFloatProperty("POSTPROCESSING_STAR_BRIGHTNESS_MIN");
            POSTPROCESSING_STAR_BRIGHTNESS_MAX = props.getFloatProperty("POSTPROCESSING_STAR_BRIGHTNESS_MAX");
            POSTPROCESSING_SPHERE_BRIGHTNESS_MIN = props.getFloatProperty("POSTPROCESSING_SPHERE_BRIGHTNESS_MIN");
            POSTPROCESSING_SPHERE_BRIGHTNESS_MAX = props.getFloatProperty("POSTPROCESSING_SPHERE_BRIGHTNESS_MAX");

            POSTPROCESSING_HUD_BRIGHTNESS_MIN = props.getFloatProperty("POSTPROCESSING_HUD_BRIGHTNESS_MIN");
            POSTPROCESSING_HUD_BRIGHTNESS_MAX = props.getFloatProperty("POSTPROCESSING_HUD_BRIGHTNESS_MAX");

            // Settings for the postprocessing shader
            POSTPROCESSING_OVERALL_BRIGHTNESS_DEF = props.getFloatProperty("POSTPROCESSING_OVERALL_BRIGHTNESS_DEF");
            POSTPROCESSING_AXES_BRIGHTNESS_DEF = props.getFloatProperty("POSTPROCESSING_AXES_BRIGHTNESS_DEF");
            POSTPROCESSING_POINT_GAS_BRIGHTNESS_DEF = props.getFloatProperty("POSTPROCESSING_POINT_GAS_BRIGHTNESS_DEF");
            POSTPROCESSING_STAR_HALO_BRIGHTNESS_DEF = props.getFloatProperty("POSTPROCESSING_STAR_HALO_BRIGHTNESS_DEF");
            POSTPROCESSING_STAR_BRIGHTNESS_DEF = props.getFloatProperty("POSTPROCESSING_STAR_BRIGHTNESS_DEF");
            POSTPROCESSING_SPHERE_BRIGHTNESS_DEF = props.getFloatProperty("POSTPROCESSING_SPHERE_BRIGHTNESS_DEF");

            GAS_COLOR_INVERTED = props.getBooleanProperty("GAS_COLOR_INVERTED");
            GAS_COLOR_BACKGROUND_INVERTED = props.getBooleanProperty("GAS_COLOR_BACKGROUND_INVERTED");
            GAS_COLOR_FROM_STARS = props.getBooleanProperty("GAS_COLOR_FROM_STARS");
            STAR_COLORS_EXAGGERATED = props.getBooleanProperty("STAR_COLORS_EXAGGERATED");

            WAITTIME_FOR_RETRY = props.getLongProperty("WAITTIME_FOR_RETRY");
            WAITTIME_FOR_MOVIE = props.getLongProperty("WAITTIME_FOR_MOVIE");

        } catch (NumberFormatException e) {
            logger.debug(e.getMessage());
        }

        currentDescription = new GlueSceneDescription(0, 0, "hotres", "");

        pointGasBlurTypeSetting = SLIDER_DEF_BLUR_TYPE;
        pointGasBlurPassSetting = SLIDER_DEF_BLUR_PASSES;
        pointGasBlurSizeSetting = SLIDER_DEF_BLUR_SIZE;

        starHaloBlurTypeSetting = SLIDER_DEF_BLUR_TYPE;
        starHaloBlurPassSetting = SLIDER_DEF_BLUR_PASSES;
        starHaloBlurSizeSetting = SLIDER_DEF_BLUR_SIZE;

    }

    public GlueSceneDescription getCurrentDescription() {
        return currentDescription;
    }

    public void setCurrentSceneDescriptionString(String value) {
        GlueSceneDescription tempDesc = currentDescription.clone();
        tempDesc.setDescriptionString(value);

        currentDescription = tempDesc;
    }

    public void setCurrentFrameNumber(int value) {
        GlueSceneDescription tempDesc = currentDescription.clone();
        tempDesc.setFrameNumber(value);

        currentDescription = tempDesc;
    }

    public void setCurrentLOD(int value) {
        GlueSceneDescription tempDesc = currentDescription.clone();
        tempDesc.setLevelOfDetail(value);

        currentDescription = tempDesc;
    }

    public void setCurrentColorMap(String value) {
        GlueSceneDescription tempDesc = currentDescription.clone();
        tempDesc.setColorMap(value);

        currentDescription = tempDesc;
    }

    public float getEpsilon() {
        return EPSILON;
    }

    public boolean getGasInvertedBackgroundColor() {
        return GAS_COLOR_BACKGROUND_INVERTED;
    }

    public boolean getGasInvertedColor() {
        return GAS_COLOR_INVERTED;
    }

    public boolean getGasStarInfluencedColor() {
        return GAS_COLOR_FROM_STARS;
    }

    public float getPostprocessingAxesBrightness() {
        return POSTPROCESSING_AXES_BRIGHTNESS_DEF;
    }

    public float getPostprocessingAxesBrightnessMax() {
        return POSTPROCESSING_AXES_BRIGHTNESS_MAX;
    }

    public float getPostprocessingAxesBrightnessMin() {
        return POSTPROCESSING_AXES_BRIGHTNESS_MIN;
    }

    public float getPostprocessingPointGasBrightness() {
        return POSTPROCESSING_POINT_GAS_BRIGHTNESS_DEF;
    }

    public float getPostprocessingPointGasBrightnessMax() {
        return POSTPROCESSING_POINT_GAS_BRIGHTNESS_MAX;
    }

    public float getPostprocessingPointGasBrightnessMin() {
        return POSTPROCESSING_POINT_GAS_BRIGHTNESS_MIN;
    }

    public float getPostprocessingHudBrightness() {
        return POSTPROCESSING_HUD_BRIGHTNESS_DEF;
    }

    public float getPostprocessingHudBrightnessMax() {
        return POSTPROCESSING_HUD_BRIGHTNESS_MAX;
    }

    public float getPostprocessingHudBrightnessMin() {
        return POSTPROCESSING_HUD_BRIGHTNESS_MIN;
    }

    public float getPostprocessingOverallBrightness() {
        return POSTPROCESSING_OVERALL_BRIGHTNESS_DEF;
    }

    public float getPostprocessingOverallBrightnessMax() {
        return POSTPROCESSING_OVERALL_BRIGHTNESS_MAX;
    }

    public float getPostprocessingOverallBrightnessMin() {
        return POSTPROCESSING_OVERALL_BRIGHTNESS_MIN;
    }

    public float getPostprocessingStarBrightness() {
        return POSTPROCESSING_STAR_BRIGHTNESS_DEF;
    }

    public float getPostprocessingStarBrightnessMax() {
        return POSTPROCESSING_STAR_BRIGHTNESS_MAX;
    }

    public float getPostprocessingStarBrightnessMin() {
        return POSTPROCESSING_STAR_BRIGHTNESS_MIN;
    }

    public float getPostprocessingStarHaloBrightness() {
        return POSTPROCESSING_STAR_HALO_BRIGHTNESS_DEF;
    }

    public float getPostprocessingStarHaloBrightnessMax() {
        return POSTPROCESSING_STAR_HALO_BRIGHTNESS_MAX;
    }

    public float getPostprocessingStarHaloBrightnessMin() {
        return POSTPROCESSING_STAR_HALO_BRIGHTNESS_MIN;
    }

    public float getPostprocessingSphereBrightnessMin() {
        return POSTPROCESSING_SPHERE_BRIGHTNESS_MIN;
    }

    public float getPostprocessingSphereBrightnessMax() {
        return POSTPROCESSING_SPHERE_BRIGHTNESS_MAX;
    }

    public float getPostprocessingSphereBrightness() {
        return POSTPROCESSING_SPHERE_BRIGHTNESS_DEF;
    }

    public float getParticleSizeMultiplierMin() {
        return PARTICLE_SIZE_MULTIPLIER_MIN;
    }

    public float getParticleSizeMultiplierMax() {
        return PARTICLE_SIZE_MULTIPLIER_MAX;
    }

    public float getParticleSizeMultiplier() {
        return PARTICLE_SIZE_MULTIPLIER_DEF;
    }

    public boolean getStarColorsExaggerated() {
        return STAR_COLORS_EXAGGERATED;
    }

    public long getWaitTimeMovie() {
        return WAITTIME_FOR_MOVIE;
    }

    public long getWaitTimeRetry() {
        return WAITTIME_FOR_RETRY;
    }

    public void setGasInvertedBackgroundColor(int stateChange) {
        if (stateChange == 1) {
            GAS_COLOR_BACKGROUND_INVERTED = true;
        }
        if (stateChange == 2) {
            GAS_COLOR_BACKGROUND_INVERTED = false;
        }
    }

    public void setInvertGasColor(int stateChange) {
        if (stateChange == 1) {
            GAS_COLOR_INVERTED = true;
        }
        if (stateChange == 2) {
            GAS_COLOR_INVERTED = false;
        }
    }

    public void setPostprocessingAxesBrightness(float value) {
        POSTPROCESSING_AXES_BRIGHTNESS_DEF = value;
    }

    public void setPostprocessingPointGasBrightness(float value) {
        POSTPROCESSING_POINT_GAS_BRIGHTNESS_DEF = value;
    }

    public void setPostprocessingHudBrightness(int value) {
        POSTPROCESSING_HUD_BRIGHTNESS_DEF = value;
    }

    public void setPostprocessingOverallBrightness(float value) {
        POSTPROCESSING_OVERALL_BRIGHTNESS_DEF = value;
    }

    public void setPostprocessingSphereBrightness(float value) {
        POSTPROCESSING_SPHERE_BRIGHTNESS_DEF = value;
    }

    public void setPostprocessingStarBrightness(float value) {
        POSTPROCESSING_STAR_BRIGHTNESS_DEF = value;
    }

    public void setPostprocessingStarHaloBrightness(float value) {
        POSTPROCESSING_STAR_HALO_BRIGHTNESS_DEF = value;
    }

    public void setParticleSizeMultiplier(float value) {
        PARTICLE_SIZE_MULTIPLIER_DEF = value;
    }

    public void setStarColorsExaggerated(int stateChange) {
        if (stateChange == 1) {
            STAR_COLORS_EXAGGERATED = true;
        }
        if (stateChange == 2) {
            STAR_COLORS_EXAGGERATED = false;
        }
    }

    public void setStarInfluencedGasColor(int stateChange) {
        if (stateChange == 1) {
            GAS_COLOR_FROM_STARS = true;
        }
        if (stateChange == 2) {
            GAS_COLOR_FROM_STARS = false;
        }
    }

    public void setWaitTimeMovie(long value) {
        WAITTIME_FOR_MOVIE = value;
    }

    public void setBezierInterpolation(boolean value) {
        BEZIER_INTERPOLATION = value;
    }

    public boolean getBezierInterpolation() {
        return BEZIER_INTERPOLATION;
    }

    public void setBezierInterpolationSteps(int value) {
        BEZIER_INTERPOLATION_STEPS = value;
    }

    public int getBezierInterpolationSteps() {
        return BEZIER_INTERPOLATION_STEPS;
    }

    public int getPreprocessAmount() {
        return PREPROCESSING_AMOUNT;
    }

    public void setPreprocessAmount(int value) {
        PREPROCESSING_AMOUNT = value;
    }

    public String getCurrentExtension() {
        return currentExtension;
    }

    public int getTimestep() {
        return timeStep;
    }

    public String[] getAcceptableExtensions() {
        return ACCEPTABLE_EXTENSIONS;
    }

    public double getStarDefaultLuminosity() {
        return STAR_DEFAULT_LUMINOSITY;
    }

    public float getMinGasDensity() {
        return 0;
    }

    public int getBlurTypeMin() {
        return SLIDER_MIN_BLUR_TYPE;
    }

    public int getBlurTypeMax() {
        return SLIDER_MAX_BLUR_TYPE;
    }

    public int getBlurPassMin() {
        return SLIDER_MIN_BLUR_PASSES;
    }

    public int getBlurPassMax() {
        return SLIDER_MAX_BLUR_PASSES;
    }

    public int getBlurSizeMin() {
        return SLIDER_MIN_BLUR_SIZE;
    }

    public int getBlurSizeMax() {
        return SLIDER_MAX_BLUR_SIZE;
    }

    @Override
    public String getScreenshotPath() {
        return SCREENSHOT_PATH;
    }

    /**
     * @return the starHaloBlurTypeSetting
     */
    public int getStarHaloBlurTypeSetting() {
        return starHaloBlurTypeSetting;
    }

    /**
     * @param starHaloBlurTypeSetting
     *            the starHaloBlurTypeSetting to set
     */
    public void setStarHaloBlurTypeSetting(int starHaloBlurTypeSetting) {
        this.starHaloBlurTypeSetting = starHaloBlurTypeSetting;
    }

    /**
     * @return the starHaloBlurPassSetting
     */
    public int getStarHaloBlurPassSetting() {
        return starHaloBlurPassSetting;
    }

    /**
     * @param starHaloBlurPassSetting
     *            the starHaloBlurPassSetting to set
     */
    public void setStarHaloBlurPassSetting(int starHaloBlurPassSetting) {
        this.starHaloBlurPassSetting = starHaloBlurPassSetting;
    }

    /**
     * @return the starHaloBlurSizeSetting
     */
    public int getStarHaloBlurSizeSetting() {
        return starHaloBlurSizeSetting;
    }

    /**
     * @param starHaloBlurSizeSetting
     *            the starHaloBlurSizeSetting to set
     */
    public void setStarHaloBlurSizeSetting(int starHaloBlurSizeSetting) {
        this.starHaloBlurSizeSetting = starHaloBlurSizeSetting;
    }

    /**
     * @param pointGasBlurTypeSetting
     *            the pointGasBlurTypeSetting to set
     */
    public void setPointGasBlurTypeSetting(int pointGasBlurTypeSetting) {
        this.pointGasBlurTypeSetting = pointGasBlurTypeSetting;
    }

    /**
     * @param pointGasBlurPassSetting
     *            the pointGasBlurPassSetting to set
     */
    public void setPointGasBlurPassSetting(int pointGasBlurPassSetting) {
        this.pointGasBlurPassSetting = pointGasBlurPassSetting;
    }

    /**
     * @return the maxOctreeDepth
     */
    @Override
    public int getMaxOctreeDepth() {
        return MAX_OCTREE_DEPTH;
    }

    /**
     * @return the pointGasBlurTypeSetting
     */
    public int getPointGasBlurTypeSetting() {
        return pointGasBlurTypeSetting;
    }

    /**
     * @return the pointGasBlurPassSetting
     */
    public int getPointGasBlurPassSetting() {
        return pointGasBlurPassSetting;
    }

    /**
     * @return the pointGasBlurSizeSetting
     */
    public int getPointGasBlurSizeSetting() {
        return pointGasBlurSizeSetting;
    }

    public void setFieldOfView(float value) {
        this.fieldOfView = value;
    }

    public float getFieldOfView() {
        return fieldOfView;
    }

    public float getZNear() {
        return zNear;
    }

    public void setZNear(float value) {
        this.zNear = value;
    }

    public float getZFar() {
        return zFar;
    }

    public void setZFar(float value) {
        this.zFar = value;
    }

    public void setPointGasPointSizeSetting(int pointGasPointSizeSetting) {
        this.pointGasPointSizeSetting = pointGasPointSizeSetting;
    }

    public float getPointGasPointSizeSetting() {
        return pointGasPointSizeSetting;
    }

    public void setPointGasBlurSizeSetting(int pointGasBlurSizeSetting) {
        this.pointGasBlurSizeSetting = pointGasBlurSizeSetting;

    }

    public float getPointGasPointSizeMin() {
        return SLIDER_MIN_POINT_SIZE;
    }

    public float getPointGasPointSizeMax() {
        return SLIDER_MAX_POINT_SIZE;
    }

    public void setPointGasPointSizeDependantOnCameraDistance(boolean pointgasSizeDependantOnCameraDistance) {
        this.pointgasSizeDependantOnCameraDistance = pointgasSizeDependantOnCameraDistance;
    }

    /**
     * Getter for pointgasSizeDependantOnCameraDistance.
     * 
     * @return the pointgasSizeDependantOnCameraDistance.
     */
    public boolean isPointgasSizeDependantOnCameraDistance() {
        return pointgasSizeDependantOnCameraDistance;
    }

    public void setDoOrbit(boolean selected) {
        this.doOrbit = selected;
    }

    public boolean isDoOrbit() {
        return doOrbit;
    }

    public void setOrbitLinkedToPlayback(boolean selected) {
        this.orbitLinkedToPlayback = selected;
    }

    public boolean isOrbitLinkedToPlayback() {
        return orbitLinkedToPlayback;
    }

    public void setOrbitSpeed(float value) {
        this.orbitSpeed = value;
    }

    public float getOrbitSpeed() {
        return orbitSpeed;
    }

    public void setXOrbitSpeed(float value) {
        this.xOrbitSpeed = value;
    }

    public void setYOrbitSpeed(float value) {
        this.yOrbitSpeed = value;
    }

    public void setZOrbitSpeed(float value) {
        this.zOrbitSpeed = value;
    }

    public float getXOrbitSpeed() {
        return xOrbitSpeed;
    }

    public float getYOrbitSpeed() {
        return yOrbitSpeed;
    }

    public float getZOrbitSpeed() {
        return zOrbitSpeed;
    }

}
