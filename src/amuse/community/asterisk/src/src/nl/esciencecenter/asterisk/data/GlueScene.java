package nl.esciencecenter.asterisk.data;

import java.nio.FloatBuffer;
import java.util.ArrayList;
import java.util.Collections;

import javax.media.opengl.GL3;

import nl.esciencecenter.asterisk.AsteriskSettings;
import nl.esciencecenter.asterisk.PointGas;
import nl.esciencecenter.asterisk.SPHGas;
import nl.esciencecenter.asterisk.Snapshot;
import nl.esciencecenter.asterisk.Sphere;
import nl.esciencecenter.asterisk.Star;
import nl.esciencecenter.asterisk.interfaces.VisualScene;
import nl.esciencecenter.asterisk.visual.DistanceComparator;
import nl.esciencecenter.asterisk.visual.PlanetModel;
import nl.esciencecenter.asterisk.visual.PointCloud;
import nl.esciencecenter.asterisk.visual.SPHOctreeNode;
import nl.esciencecenter.asterisk.visual.SphereModel;
import nl.esciencecenter.asterisk.visual.StarModel;
import nl.esciencecenter.esight.exceptions.InverseNotAvailableException;
import nl.esciencecenter.esight.exceptions.UninitializedException;
import nl.esciencecenter.esight.math.MatF4;
import nl.esciencecenter.esight.math.MatrixFMath;
import nl.esciencecenter.esight.math.VecF3;
import nl.esciencecenter.esight.math.VecF4;
import nl.esciencecenter.esight.models.Model;
import nl.esciencecenter.esight.shaders.ShaderProgram;

import com.jogamp.common.nio.Buffers;

public class GlueScene implements Runnable, VisualScene {
    private final GlueSceneStorage sceneStore;
    private final GlueSceneDescription description;
    private final Snapshot scene;

    private ArrayList<StarModel> stars;
    private ArrayList<PlanetModel> planets;
    private ArrayList<SphereModel> spheres;

    private PointCloud gasParticles;
    private SPHOctreeNode gasOctree;

    private boolean initialized = false;

    public GlueScene(GlueSceneStorage sceneStore,
            GlueSceneDescription description, Snapshot scene) {
        this.sceneStore = sceneStore;
        this.description = description;
        this.scene = scene;
    }

    @Override
    public void run() {

        Star[] glueStars = scene.getStars();
        Model starBaseModel = sceneStore.getStarBaseModel();
        if (glueStars != null) {
            this.stars = new ArrayList<StarModel>();
            for (Star glueStar : glueStars) {
                StarModel starModel = new StarModel(starBaseModel, glueStar);
                stars.add(starModel);
            }
        }

        // Planet[] gluePlanets = scene.getPlanets();
        // Model planetBaseModel = sceneStore.getPlanetBaseModel();
        // if (gluePlanets != null) {
        // this.planets = new ArrayList<PlanetModel>();
        // for (Planet gluePlanet : gluePlanets) {
        // PlanetModel planetModel = new PlanetModel(planetBaseModel,
        // gluePlanet);
        // planets.add(planetModel);
        // }
        // }

        Sphere[] glueSpheres = scene.getSpheres();
        Model sphereBaseModel = sceneStore.getSphereBaseModel();
        if (glueSpheres != null) {
            this.spheres = new ArrayList<SphereModel>();
            for (Sphere glueSphere : glueSpheres) {
                SphereModel sphereModel = new SphereModel(sphereBaseModel,
                        glueSphere);
                spheres.add(sphereModel);
            }
        }

        PointGas[] pointGasses = scene.getPointGas();
        if (pointGasses != null) {
            int numPointGasParticles = pointGasses.length;

            FloatBuffer pointGasCoords = Buffers
                    .newDirectFloatBuffer(numPointGasParticles * 3);
            FloatBuffer pointGasColors = Buffers
                    .newDirectFloatBuffer(numPointGasParticles * 4);

            for (PointGas gluePointGas : pointGasses) {
                float[] coords = gluePointGas.getCoordinates();
                float[] color = gluePointGas.getColor();

                for (int i = 0; i < 3; i++) {
                    pointGasCoords.put(coords[i]);
                }

                for (int i = 0; i < 4; i++) {
                    pointGasColors.put(color[i]);
                }
            }

            gasParticles = new PointCloud(numPointGasParticles, pointGasCoords,
                    pointGasColors);
        }

        // PointGas[] pointGasses = scene.getPointGas();
        // if (pointGasses != null) {
        // VectorList coords = new VectorList(3, Vector.Type.FLOAT);
        // VectorList colors = new VectorList(4, Vector.Type.FLOAT);
        //
        // for (PointGas gluePointGas : pointGasses) {
        // float[] v = gluePointGas.getCoordinates();
        // coords.add(new VecF3(v[0], v[1], v[2]));
        //
        // float[] c = gluePointGas.getColor();
        // colors.add(new VecF4(c[0], c[1], c[2], c[3]));
        // }
        //
        // gasParticles = new PointCloud(coords, colors);
        // }

        SPHGas[] sphGasses = scene.getSphGas();
        if (sphGasses != null) {
            Model octreeBaseModel = sceneStore.getSPHOctreeBaseModel();
            this.gasOctree = new SPHOctreeNode(octreeBaseModel, 0, new VecF3(),
                    AsteriskSettings.INITIAL_OCTREE_SIZE, description);
            for (SPHGas glueSPHGas : sphGasses) {
                float[] rawCoords = glueSPHGas.getCoordinates();
                float[] rawColor = glueSPHGas.getColor();

                VecF3 coords = new VecF3(rawCoords[0], rawCoords[1],
                        rawCoords[2]);
                VecF4 color = new VecF4(rawColor[0], rawColor[1], rawColor[2],
                        rawColor[3]);

                this.gasOctree.addElement(coords, color);
            }
        }

        sceneStore.setScene(description, this);
    }

    @Override
    public synchronized void drawStars(GL3 gl, ShaderProgram program,
            MatF4 MVMatrix) {
        MatF4 viewModel;
        if (stars != null) {
            try {
                viewModel = MatrixFMath.inverse(MVMatrix);
                VecF3 cameraPos = new VecF3(viewModel.get(3), viewModel.get(7),
                        viewModel.get(11));
                Collections.sort(stars, new DistanceComparator(cameraPos));
            } catch (InverseNotAvailableException e) {
                e.printStackTrace();
            }

            for (StarModel s : stars) {
                s.draw(gl, program, MVMatrix);
            }
        }
    }

    public synchronized void drawPlanets(GL3 gl, ShaderProgram program,
            MatF4 MVMatrix) {
        if (planets != null) {
            for (PlanetModel p : planets) {
                p.draw(gl, program, MVMatrix);
            }
        }
    }

    @Override
    public synchronized void drawSpheres(GL3 gl, ShaderProgram program,
            MatF4 MVMatrix) {
        if (spheres != null) {
            for (SphereModel s : spheres) {
                s.draw(gl, program, MVMatrix);
            }
        }
    }

    @Override
    public synchronized void drawGasPointCloud(GL3 gl, ShaderProgram program,
            MatF4 MVMatrix) {
        VecF3 cameraPos;
        try {
            MatF4 viewModel = MatrixFMath.inverse(MVMatrix);
            cameraPos = new VecF3(viewModel.get(3), viewModel.get(7),
                    viewModel.get(11));
        } catch (InverseNotAvailableException e) {
            cameraPos = new VecF3();
            e.printStackTrace();
        }

        if (gasParticles != null) {
            program.setUniformVector("CameraPos", cameraPos);
            program.setUniformMatrix("MVMatrix", MVMatrix);

            try {
                program.use(gl);
            } catch (UninitializedException e) {
                e.printStackTrace();
            }

            gasParticles.draw(gl, program);
        }
    }

    @Override
    public synchronized void drawGasOctree(GL3 gl, ShaderProgram program,
            MatF4 MVMatrix) {
        if (gasOctree != null) {
            program.setUniformMatrix("MVMatrix", MVMatrix);

            program.setUniformMatrix("NormalMatrix",
                    MatrixFMath.getNormalMatrix(MVMatrix));

            gasOctree.draw(gl, program);
        }
    }

    public void init(GL3 gl) {
        if (!initialized) {
            if (stars != null) {
                for (StarModel s : stars) {
                    s.init();
                }
            }

            if (planets != null) {
                for (PlanetModel p : planets) {
                    p.init();
                }
            }

            if (spheres != null) {
                for (SphereModel s : spheres) {
                    s.init();
                }
            }

            if (gasParticles != null) {
                gasParticles.init(gl);
            }

            if (gasOctree != null) {
                gasOctree.init();
            }

            initialized = true;
        }
    }

    @Override
    public void dispose(GL3 gl) {
        gasParticles.dispose(gl);

        initialized = false;
    }

    @Override
    public String getDescriptionString() {
        return scene.getSceneDecriptionString();
    }
}
