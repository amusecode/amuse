package nl.esciencecenter.asterisk;

import nl.esciencecenter.asterisk.data.GlueSceneDescription;
import nl.esciencecenter.esight.util.Settings;
import nl.esciencecenter.esight.util.TypedProperties;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class AsteriskSettings extends Settings {
    private final static Logger logger = LoggerFactory
            .getLogger(AsteriskSettings.class);

    private final double STAR_DEFAULT_LUMINOSITY = 3000.0;

    private static class SingletonHolder {
        public final static AsteriskSettings instance = new AsteriskSettings();
    }

    private final String[] ACCEPTABLE_EXTENSIONS = new String[] { "gas", "bin" };

    private final float MIN_GAS_DENSITY = 0f;
    private final float MAX_GAS_DENSITY = 15f;

    // Minimum and maximum values for the brightness sliders
    private float POSTPROCESSING_OVERALL_BRIGHTNESS_MIN = 0f;
    private float POSTPROCESSING_OVERALL_BRIGHTNESS_MAX = 10f;
    private float POSTPROCESSING_AXES_BRIGHTNESS_MIN = 0f;
    private float POSTPROCESSING_AXES_BRIGHTNESS_MAX = 10f;
    private float POSTPROCESSING_POINT_GAS_BRIGHTNESS_MIN = 0f;
    private float POSTPROCESSING_POINT_GAS_BRIGHTNESS_MAX = 10f;
    private float POSTPROCESSING_OCTREE_GAS_BRIGHTNESS_MIN = 0f;
    private float POSTPROCESSING_OCTREE_GAS_BRIGHTNESS_MAX = 10f;
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
    private float POSTPROCESSING_POINT_GAS_BRIGHTNESS_DEF = 4f;
    private float POSTPROCESSING_OCTREE_GAS_BRIGHTNESS_DEF = 4f;
    private float POSTPROCESSING_STAR_HALO_BRIGHTNESS_DEF = 4f;
    private float POSTPROCESSING_STAR_BRIGHTNESS_DEF = 4f;
    private float POSTPROCESSING_SPHERE_BRIGHTNESS_DEF = 4f;
    private float POSTPROCESSING_HUD_BRIGHTNESS_DEF = 4f;

    private final float PARTICLE_SIZE_MULTIPLIER_MIN = 1f;
    private float PARTICLE_SIZE_MULTIPLIER_DEF = 10f;
    private final float PARTICLE_SIZE_MULTIPLIER_MAX = 100f;

    private final int SLIDER_MIN_BLUR_TYPE = 1;
    private final int SLIDER_DEF_BLUR_TYPE = 1;
    private final int SLIDER_MAX_BLUR_TYPE = 8;

    private final int SLIDER_MIN_BLUR_PASSES = 0;
    private final int SLIDER_DEF_BLUR_PASSES = 1;
    private final int SLIDER_MAX_BLUR_PASSES = 2;

    private final int SLIDER_MIN_BLUR_SIZE = 2;
    private final int SLIDER_DEF_BLUR_SIZE = 2;
    private final int SLIDER_MAX_BLUR_SIZE = 8;

    private boolean GAS_COLOR_INVERTED = false;
    private boolean GAS_COLOR_BACKGROUND_INVERTED = false;
    private boolean GAS_COLOR_FROM_STARS = false;

    private boolean STAR_COLORS_EXAGGERATED = true;

    private long WAITTIME_FOR_RETRY = 10000;
    private long WAITTIME_FOR_MOVIE = 200;
    private float EPSILON = 1.0E-7f;

    private float GAS_OPACITY_FACTOR_MIN = 0f;
    private float GAS_OPACITY_FACTOR_DEF = 1f;
    private float GAS_OPACITY_FACTOR_MAX = 2f;

    private static final int GAS_OPACITY_RATIO_MIN = 1;
    private static final int GAS_OPACITY_RATIO_MAX = 100;

    private static final boolean OCTREE_RANDOM_OFFSET = false;

    private boolean BEZIER_INTERPOLATION = false;
    private int BEZIER_INTERPOLATION_STEPS = 10;

    private int PREPROCESSING_AMOUNT = 5;

    private int pointGasBlurTypeSetting;
    private int pointGasBlurPassSetting;
    private int pointGasBlurSizeSetting;

    private int octreeGasBlurTypeSetting;
    private int octreeGasBlurPassSetting;
    private int octreeGasBlurSizeSetting;

    private int starHaloBlurTypeSetting;
    private int starHaloBlurPassSetting;
    private int starHaloBlurSizeSetting;

    private float gasOpacityRatio = 10;

    public static int MAX_OCTREE_ELEMENTS_PER_NODE = 2;
    public static final int MAX_OCTREE_DEPTH = 25;

    public static final float INITIAL_OCTREE_SIZE = 500f;

    public static final int STAR_SUBDIVISION = 2;
    public static final int PLANET_SUBDIVISION = 2;
    public static final int SPHERE_SUBDIVISION = 2;

    public static final int OCTREE_MODEL_SUBDIVISION = 2;
    
    private float fieldOfView = 45.0f;
    private float zNear = 0.1f;
    private float zFar = 3000.0f;

    private GlueSceneDescription currentDescription;

    public static AsteriskSettings getInstance() {
        return SingletonHolder.instance;
    }

    private final String currentExtension = "bin";
    private final int timeStep = 1;

    private AsteriskSettings() {
        super();

        try {
            final TypedProperties props = new TypedProperties();
            props.loadFromClassPath("settings.properties");

            // Minimum and maximum values for the brightness sliders
            POSTPROCESSING_OVERALL_BRIGHTNESS_MIN = props
                    .getFloatProperty("POSTPROCESSING_OVERALL_BRIGHTNESS_MIN");
            POSTPROCESSING_OVERALL_BRIGHTNESS_MAX = props
                    .getFloatProperty("POSTPROCESSING_OVERALL_BRIGHTNESS_MAX");
            POSTPROCESSING_AXES_BRIGHTNESS_MIN = props
                    .getFloatProperty("POSTPROCESSING_AXES_BRIGHTNESS_MIN");
            POSTPROCESSING_AXES_BRIGHTNESS_MAX = props
                    .getFloatProperty("POSTPROCESSING_AXES_BRIGHTNESS_MAX");
            POSTPROCESSING_POINT_GAS_BRIGHTNESS_MIN = props
                    .getFloatProperty("POSTPROCESSING_POINT_GAS_BRIGHTNESS_MIN");
            POSTPROCESSING_POINT_GAS_BRIGHTNESS_MAX = props
                    .getFloatProperty("POSTPROCESSING_POINT_GAS_BRIGHTNESS_MAX");
            POSTPROCESSING_OCTREE_GAS_BRIGHTNESS_MIN = props
                    .getFloatProperty("POSTPROCESSING_OCTREE_GAS_BRIGHTNESS_MIN");
            POSTPROCESSING_OCTREE_GAS_BRIGHTNESS_MAX = props
                    .getFloatProperty("POSTPROCESSING_OCTREE_GAS_BRIGHTNESS_MAX");
            POSTPROCESSING_STAR_HALO_BRIGHTNESS_MIN = props
                    .getFloatProperty("POSTPROCESSING_STAR_HALO_BRIGHTNESS_MIN");
            POSTPROCESSING_STAR_HALO_BRIGHTNESS_MAX = props
                    .getFloatProperty("POSTPROCESSING_STAR_HALO_BRIGHTNESS_MAX");
            POSTPROCESSING_STAR_BRIGHTNESS_MIN = props
                    .getFloatProperty("POSTPROCESSING_STAR_BRIGHTNESS_MIN");
            POSTPROCESSING_STAR_BRIGHTNESS_MAX = props
                    .getFloatProperty("POSTPROCESSING_STAR_BRIGHTNESS_MAX");
            POSTPROCESSING_SPHERE_BRIGHTNESS_MIN = props
                    .getFloatProperty("POSTPROCESSING_SPHERE_BRIGHTNESS_MIN");
            POSTPROCESSING_SPHERE_BRIGHTNESS_MAX = props
                    .getFloatProperty("POSTPROCESSING_SPHERE_BRIGHTNESS_MAX");

            POSTPROCESSING_HUD_BRIGHTNESS_MIN = props
                    .getFloatProperty("POSTPROCESSING_HUD_BRIGHTNESS_MIN");
            POSTPROCESSING_HUD_BRIGHTNESS_MAX = props
                    .getFloatProperty("POSTPROCESSING_HUD_BRIGHTNESS_MAX");

            // Settings for the postprocessing shader
            POSTPROCESSING_OVERALL_BRIGHTNESS_DEF = props
                    .getFloatProperty("POSTPROCESSING_OVERALL_BRIGHTNESS_DEF");
            POSTPROCESSING_AXES_BRIGHTNESS_DEF = props
                    .getFloatProperty("POSTPROCESSING_AXES_BRIGHTNESS_DEF");
            POSTPROCESSING_POINT_GAS_BRIGHTNESS_DEF = props
                    .getFloatProperty("POSTPROCESSING_POINT_GAS_BRIGHTNESS_DEF");
            POSTPROCESSING_OCTREE_GAS_BRIGHTNESS_DEF = props
                    .getFloatProperty("POSTPROCESSING_OCTREE_GAS_BRIGHTNESS_DEF");
            POSTPROCESSING_STAR_HALO_BRIGHTNESS_DEF = props
                    .getFloatProperty("POSTPROCESSING_STAR_HALO_BRIGHTNESS_DEF");
            POSTPROCESSING_STAR_BRIGHTNESS_DEF = props
                    .getFloatProperty("POSTPROCESSING_STAR_BRIGHTNESS_DEF");
            POSTPROCESSING_SPHERE_BRIGHTNESS_DEF = props
                    .getFloatProperty("POSTPROCESSING_SPHERE_BRIGHTNESS_DEF");

            GAS_COLOR_INVERTED = props.getBooleanProperty("GAS_COLOR_INVERTED");
            GAS_COLOR_BACKGROUND_INVERTED = props
                    .getBooleanProperty("GAS_COLOR_BACKGROUND_INVERTED");
            GAS_COLOR_FROM_STARS = props
                    .getBooleanProperty("GAS_COLOR_FROM_STARS");
            STAR_COLORS_EXAGGERATED = props
                    .getBooleanProperty("STAR_COLORS_EXAGGERATED");

            WAITTIME_FOR_RETRY = props.getLongProperty("WAITTIME_FOR_RETRY");
            WAITTIME_FOR_MOVIE = props.getLongProperty("WAITTIME_FOR_MOVIE");
            EPSILON = props.getFloatProperty("EPSILON");

            GAS_OPACITY_FACTOR_MIN = props
                    .getFloatProperty("GAS_OPACITY_FACTOR_MIN");
            GAS_OPACITY_FACTOR_DEF = props
                    .getFloatProperty("GAS_OPACITY_FACTOR_DEF");
            GAS_OPACITY_FACTOR_MAX = props
                    .getFloatProperty("GAS_OPACITY_FACTOR_MAX");

            BEZIER_INTERPOLATION = props
                    .getBooleanProperty("BEZIER_INTERPOLATION");
            BEZIER_INTERPOLATION_STEPS = props
                    .getIntProperty("BEZIER_INTERPOLATION_STEPS");
            MAX_OCTREE_ELEMENTS_PER_NODE = props
                    .getIntProperty("MAX_OCTREE_ELEMENTS_PER_NODE");

            PREPROCESSING_AMOUNT = props.getIntProperty("PREPROCESSING_AMOUNT");

        } catch (NumberFormatException e) {
            logger.debug(e.getMessage());
        }

        currentDescription = new GlueSceneDescription(0, 0, "hotres",
                MIN_GAS_DENSITY, MAX_GAS_DENSITY, "");

        // currentDescription = new GlueSceneDescription(0, 0, "default", 0f,
        // 25f);

        pointGasBlurTypeSetting = SLIDER_DEF_BLUR_TYPE;
        pointGasBlurPassSetting = SLIDER_DEF_BLUR_PASSES;
        pointGasBlurSizeSetting = SLIDER_DEF_BLUR_SIZE;

        octreeGasBlurTypeSetting = SLIDER_DEF_BLUR_TYPE;
        octreeGasBlurPassSetting = SLIDER_DEF_BLUR_PASSES;
        octreeGasBlurSizeSetting = SLIDER_DEF_BLUR_SIZE;

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

    public void setCurrentLowerBound(float value) {
        GlueSceneDescription tempDesc = currentDescription.clone();
        tempDesc.setLowerBound(value);

        currentDescription = tempDesc;
    }

    public void setCurrentUpperBound(float value) {
        GlueSceneDescription tempDesc = currentDescription.clone();
        tempDesc.setUpperBound(value);

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

    public float getPostprocessingOctreeGasBrightness() {
        return POSTPROCESSING_OCTREE_GAS_BRIGHTNESS_DEF;
    }

    public float getPostprocessingOctreeGasBrightnessMax() {
        return POSTPROCESSING_OCTREE_GAS_BRIGHTNESS_MAX;
    }

    public float getPostprocessingOctreeGasBrightnessMin() {
        return POSTPROCESSING_OCTREE_GAS_BRIGHTNESS_MIN;
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

    public void setPostprocessingOctreeGasBrightness(float value) {
        POSTPROCESSING_OCTREE_GAS_BRIGHTNESS_DEF = value;
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

    public void setGasOpacityFactor(float value) {
        GAS_OPACITY_FACTOR_DEF = value;
    }

    public float getGasOpacityFactorMin() {
        return GAS_OPACITY_FACTOR_MIN;
    }

    public float getGasOpacityFactor() {
        return GAS_OPACITY_FACTOR_DEF;
    }

    public float getGasOpacityFactorMax() {
        return GAS_OPACITY_FACTOR_MAX;
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

    public void setVariableRange(int sliderLowerValue, int sliderUpperValue) {
        float diff = MAX_GAS_DENSITY - MIN_GAS_DENSITY;
        float minFloatValue = (sliderLowerValue / 100f) * diff;
        float maxFloatValue = (sliderUpperValue / 100f) * diff;

        GlueSceneDescription result = currentDescription.clone();
        result.setLowerBound(minFloatValue);
        result.setUpperBound(maxFloatValue);
        currentDescription = result;
    }

    public int getRangeSliderLowerValue() {
        float min = MIN_GAS_DENSITY;
        float max = MAX_GAS_DENSITY;
        float currentMin = currentDescription.getLowerBound();

        float diff = max - min;
        float result = (currentMin - min) / diff;

        return (int) (result * 100) - 1;
    }

    public int getRangeSliderUpperValue() {
        float min = MIN_GAS_DENSITY;
        float max = MAX_GAS_DENSITY;
        float currentMax = currentDescription.getUpperBound();

        float diff = max - min;
        float result = (currentMax - min) / diff;

        return (int) (result * 100) - 1;
    }

    public float getMinGasDensity() {
        return 0;
    }

    public float getMaxGasDensity() {
        return MAX_GAS_DENSITY;
    }

    public float getGasOpacityRatio() {
        return gasOpacityRatio;
    }

    public void setGasOpacityRatio(float value) {
        gasOpacityRatio = value;
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

    public int getGasOpacityRatioMin() {
        return GAS_OPACITY_RATIO_MIN;
    }

    public int getGasOpacityRatioMax() {
        return GAS_OPACITY_RATIO_MAX;
    }

    public boolean isGasRandomCenterOffset() {
        return OCTREE_RANDOM_OFFSET;
    }

    @Override
    public String getScreenshotPath() {
        return SCREENSHOT_PATH;
    }

    /**
     * @return the octreeGasBlurTypeSetting
     */
    public int getOctreeGasBlurTypeSetting() {
        return octreeGasBlurTypeSetting;
    }

    /**
     * @param octreeGasBlurTypeSetting
     *            the octreeGasBlurTypeSetting to set
     */
    public void setOctreeGasBlurTypeSetting(int octreeGasBlurTypeSetting) {
        this.octreeGasBlurTypeSetting = octreeGasBlurTypeSetting;
    }

    /**
     * @return the octreeGasBlurPassSetting
     */
    public int getOctreeGasBlurPassSetting() {
        return octreeGasBlurPassSetting;
    }

    /**
     * @param octreeGasBlurPassSetting
     *            the octreeGasBlurPassSetting to set
     */
    public void setOctreeGasBlurPassSetting(int octreeGasBlurPassSetting) {
        this.octreeGasBlurPassSetting = octreeGasBlurPassSetting;
    }

    /**
     * @return the octreeGasBlurSizeSetting
     */
    public int getOctreeGasBlurSizeSetting() {
        return octreeGasBlurSizeSetting;
    }

    /**
     * @param octreeGasBlurSizeSetting
     *            the octreeGasBlurSizeSetting to set
     */
    public void setOctreeGasBlurSizeSetting(int octreeGasBlurSizeSetting) {
        this.octreeGasBlurSizeSetting = octreeGasBlurSizeSetting;
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
     * @param pointGasBlurSizeSetting
     *            the pointGasBlurSizeSetting to set
     */
    public void setPointGasBlurSizeSetting(int pointGasBlurSizeSetting) {
        this.pointGasBlurSizeSetting = pointGasBlurSizeSetting;
    }

    /**
     * @return the maxOctreeElementsPerNode
     */
    public static int getMaxOctreeElementsPerNode() {
        return MAX_OCTREE_ELEMENTS_PER_NODE;
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

    
}
