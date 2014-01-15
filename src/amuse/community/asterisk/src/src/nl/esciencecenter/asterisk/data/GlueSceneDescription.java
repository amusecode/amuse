package nl.esciencecenter.asterisk.data;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class GlueSceneDescription {
    private final static Logger logger = LoggerFactory.getLogger(GlueSceneDescription.class);

    private int frameNumber;
    private String colorMap;
    private int levelOfDetail;

    private String descriptionString;

    public GlueSceneDescription(int frameNumber, int levelOfDetail, String colorMap, String sceneDescriptionString) {
        this.frameNumber = frameNumber;
        this.levelOfDetail = levelOfDetail;
        this.colorMap = colorMap;
        this.descriptionString = sceneDescriptionString;
    }

    @Override
    public int hashCode() {
        int frameNumberPrime = (frameNumber + 131) * 1543;
        int lodPrime = (frameNumber + 439) * 1847;
        int colorMapPrime = (colorMap.hashCode() + 919) * 7883;
        int descriptionPrime = descriptionString.hashCode();

        int hashCode = frameNumberPrime + lodPrime + colorMapPrime + descriptionPrime;

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
        return (frameNumber == that.frameNumber && levelOfDetail == that.levelOfDetail
                && colorMap.compareTo(that.colorMap) == 0 && descriptionString.compareTo(that.descriptionString) == 0);
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

    @Override
    public GlueSceneDescription clone() {
        return new GlueSceneDescription(frameNumber, levelOfDetail, colorMap, descriptionString);
    }

    public void setFrameNumber(int frameNumber) {
        this.frameNumber = frameNumber;
    }

    public void setColorMap(String colorMap) {
        this.colorMap = colorMap;
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
