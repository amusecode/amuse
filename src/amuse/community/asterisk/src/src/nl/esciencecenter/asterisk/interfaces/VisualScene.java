package nl.esciencecenter.asterisk.interfaces;

import javax.media.opengl.GL3;

import nl.esciencecenter.esight.math.MatF4;
import nl.esciencecenter.esight.shaders.ShaderProgram;

public interface VisualScene {

    public void dispose(GL3 gl);

    public void drawGasPointCloud(GL3 gl, ShaderProgram program, MatF4 mv);

    public void drawSpheres(GL3 gl, ShaderProgram program, MatF4 mv);

    public void drawStars(GL3 gl, ShaderProgram program, MatF4 mv);

    public void drawGasOctree(GL3 gl, ShaderProgram program, MatF4 mv);

    public String getDescriptionString();

}
