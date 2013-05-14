package nl.esciencecenter.asterisk.visual;

import java.util.Random;

import javax.media.opengl.GL3;

import nl.esciencecenter.asterisk.Star;
import nl.esciencecenter.esight.exceptions.InverseNotAvailableException;
import nl.esciencecenter.esight.exceptions.UninitializedException;
import nl.esciencecenter.esight.math.MatF4;
import nl.esciencecenter.esight.math.MatrixFMath;
import nl.esciencecenter.esight.math.VecF3;
import nl.esciencecenter.esight.math.VecF4;
import nl.esciencecenter.esight.models.Model;
import nl.esciencecenter.esight.shaders.ShaderProgram;

public class StarModel {
    private final Model baseModel;
    private final Star glueStar;

    private VecF3 coords;
    private VecF4 color;
    private float radius;
    private float offsetRandomValue;

    private boolean initialized;

    public StarModel(Model baseModel, Star glueStar) {
        this.baseModel = baseModel;
        this.glueStar = glueStar;

        initialized = false;
    }

    public void init() {
        if (!initialized) {
            float[] rawCoords = glueStar.getCoordinates();
            float[] rawColor = glueStar.getColor();
            int index = glueStar.getIndex();

            offsetRandomValue = (new Random(index * 10000)).nextFloat() * 3.14f;

            // System.err.println("Offset Random Value: " + offsetRandomValue);

            coords = new VecF3(rawCoords[0], rawCoords[1], rawCoords[2]);
            color = new VecF4(rawColor[0], rawColor[1], rawColor[2],
                    rawColor[3]);
            radius = glueStar.getRadius();

            initialized = true;
        }
    }

    public void draw(GL3 gl, ShaderProgram program, MatF4 MVMatrix) {
        if (!initialized) {
            init();
        }

        program.setUniformMatrix("ScaleMatrix", MatrixFMath.scale(radius));

        MatF4 newM = MVMatrix.mul(MatrixFMath.translate(coords));
        program.setUniformMatrix("MVMatrix", newM);

        program.setUniformMatrix("NormalMatrix",
                MatrixFMath.getNormalMatrix(newM));

        VecF3 cameraPos;
        try {
            MatF4 viewModel = MatrixFMath.inverse(newM);
            cameraPos = new VecF3(viewModel.get(3), viewModel.get(7),
                    viewModel.get(11));
        } catch (InverseNotAvailableException e) {
            cameraPos = new VecF3();
            e.printStackTrace();
        }

        program.setUniformVector("CameraPos", cameraPos);

        program.setUniformVector("Color", color);
        program.setUniform("OffsetRandomValue", offsetRandomValue);

        try {
            program.use(gl);
        } catch (UninitializedException e) {
            e.printStackTrace();
        }
        baseModel.draw(gl, program);
    }

    public VecF4 getColor() {
        if (!initialized) {
            init();
        }
        return color;
    }

    public float getRadius() {
        if (!initialized) {
            init();
        }

        return radius;
    }

    public VecF3 getLocation() {
        if (!initialized) {
            init();
        }
        return coords;
    }
}
