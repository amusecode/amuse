package nl.esciencecenter.asterisk.data;

import java.util.HashMap;

import javax.media.opengl.GL3;

import nl.esciencecenter.asterisk.AsteriskSettings;
import nl.esciencecenter.asterisk.interfaces.SceneStorage;
import nl.esciencecenter.asterisk.interfaces.VisualScene;
import nl.esciencecenter.esight.models.Model;
import nl.esciencecenter.esight.models.Sphere;

public class GlueSceneStorage implements SceneStorage {
    private GlueSceneDescription oldScene = null;
    private GlueSceneDescription newScene = null;
    private HashMap<GlueSceneDescription, GlueScene> sceneStorage;

    private final GlueDatasetManager manager;

    private final Sphere starBaseModel;
    private final Sphere planetBaseModel;
    private final Sphere sphereBaseModel;

    private boolean initialized = false;

    public GlueSceneStorage(GlueDatasetManager manager) {
        sceneStorage = new HashMap<GlueSceneDescription, GlueScene>();

        this.starBaseModel = new Sphere(AsteriskSettings.STAR_SUBDIVISION, true);
        this.planetBaseModel = new Sphere(AsteriskSettings.PLANET_SUBDIVISION, false);
        this.sphereBaseModel = new Sphere(AsteriskSettings.SPHERE_SUBDIVISION, false);

        this.manager = manager;
    }

    @Override
    public void init(GL3 gl) {
        if (!initialized) {
            starBaseModel.init(gl);
            planetBaseModel.init(gl);
            sphereBaseModel.init(gl);

            initialized = true;
        }
    }

    @Override
    public synchronized VisualScene getScene() {
        GlueScene result = null;

        if (sceneStorage.containsKey(newScene)) {
            sceneStorage.remove(oldScene);

            result = sceneStorage.get(newScene);
        } else {
            result = sceneStorage.get(oldScene);
        }

        return result;
    }

    @Override
    public synchronized void requestNewConfiguration(GlueSceneDescription newDescription) {
        // System.out.println("New config request");
        HashMap<GlueSceneDescription, GlueScene> newSceneStore = new HashMap<GlueSceneDescription, GlueScene>();

        for (GlueSceneDescription description : sceneStorage.keySet()) {
            if (description == oldScene || description == newScene) {
                newSceneStore.put(description, sceneStorage.get(description));
            }
        }
        sceneStorage = newSceneStore;

        oldScene = newScene;
        newScene = newDescription;
        if (!sceneStorage.containsValue(newScene)) {
            manager.buildScene(newScene);
        }
    }

    public boolean doneWithLastRequest() {
        boolean failure = false;

        if (sceneStorage.get(newScene) == null) {
            failure = true;
        }

        return !failure;
    }

    public void setScene(GlueSceneDescription description, GlueScene scene) {
        sceneStorage.put(description, scene);
    }

    public Model getStarBaseModel() {
        return starBaseModel;
    }

    public Model getPlanetBaseModel() {
        return planetBaseModel;
    }

    public Model getSphereBaseModel() {
        return sphereBaseModel;
    }
}
