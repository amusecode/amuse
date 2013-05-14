package nl.esciencecenter.asterisk.interfaces;

import javax.media.opengl.GL3;

import nl.esciencecenter.asterisk.data.GlueSceneDescription;

public interface SceneStorage {

    public void init(GL3 gl);

    public void requestNewConfiguration(GlueSceneDescription currentDescription);

    public VisualScene getScene();

}
