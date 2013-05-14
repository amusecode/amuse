package nl.esciencecenter.asterisk.data;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class GlueSceneDescription {
    private final static Logger logger = LoggerFactory
            .getLogger(GlueSceneDescription.class);

    private int frameNumber;
    private String colorMap;
    private float lowerBound;
    private float upperBound;
    private int levelOfDetail;

    private String descriptionString;

    public GlueSceneDescription(int frameNumber, int levelOfDetail,
            String colorMap, float lowerBound, float upperBound,
            String sceneDescriptionString) {
        this.frameNumber = frameNumber;
        this.levelOfDetail = levelOfDetail;
        this.colorMap = colorMap;
        this.lowerBound = lowerBound;
        this.upperBound = upperBound;
        this.descriptionString = sceneDescriptionString;
    }

    @Override
    public int hashCode() {
        int frameNumberPrime = (frameNumber + 131) * 1543;
        int lodPrime = (frameNumber + 439) * 1847;
        int colorMapPrime = (colorMap.hashCode() + 919) * 7883;
        int lowerBoundPrime = (int) ((lowerBound + 41) * 1543);
        int upperBoundPrime = (int) ((upperBound + 67) * 2957);
        int descriptionPrime = descriptionString.hashCode();

        int hashCode = frameNumberPrime + lodPrime + colorMapPrime
                + lowerBoundPrime + upperBoundPrime + descriptionPrime;

        return hashCode;
    }

    @Override
    public boolean equals(Object thatObject) {
        if (this == thatObject)
            return true;
        if (!(thatObject instanceof GlueSceneDescription))
            return false;

        // cast to native object is now safe
        GlueSceneDescription that = (GlueSceneDescription) thatObject;

        // now a proper field-by-field evaluation can be made
        return (frameNumber == that.frameNumber
                && levelOfDetail == that.levelOfDetail
                && lowerBound == that.lowerBound
                && upperBound == that.upperBound
                && colorMap.compareTo(that.colorMap) == 0 && descriptionString
                    .compareTo(that.descriptionString) == 0);
    }

    public static Logger getLogger() {
        return logger;
    }

    public int getFrameNumber() {
        return frameNumber;
    }

    public String getColorMap() {
        return colorMap;
    }

    public float getLowerBound() {
        return lowerBound;
    }

    public float getUpperBound() {
        return upperBound;
    }

    @Override
    public GlueSceneDescription clone() {
        return new GlueSceneDescription(frameNumber, levelOfDetail, colorMap,
                lowerBound, upperBound, descriptionString);
    }

    public void setFrameNumber(int frameNumber) {
        this.frameNumber = frameNumber;
    }

    public void setColorMap(String colorMap) {
        this.colorMap = colorMap;
    }

    public void setLowerBound(float lowerBound) {
        this.lowerBound = lowerBound;
    }

    public void setUpperBound(float upperBound) {
        this.upperBound = upperBound;
    }

    public int getLevelOfDetail() {
        return levelOfDetail;
    }

    public void setLevelOfDetail(int value) {
        levelOfDetail = value;
    }

    public String getDescriptionString() {
        return descriptionString;
    }

    public void setDescriptionString(String descriptionString) {
        this.descriptionString = descriptionString;
    }
}
