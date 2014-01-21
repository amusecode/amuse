package nl.esciencecenter.asterisk;

import java.io.File;
import java.util.concurrent.TimeUnit;

import javax.media.opengl.GL;
import javax.media.opengl.GL2GL3;
import javax.media.opengl.GL3;
import javax.media.opengl.GLAutoDrawable;
import javax.media.opengl.GLContext;
import javax.media.opengl.GLEventListener;
import javax.media.opengl.GLException;

import nl.esciencecenter.asterisk.data.GlueScene;
import nl.esciencecenter.asterisk.data.GlueSceneDescription;
import nl.esciencecenter.asterisk.data.GlueTimedPlayer;
import nl.esciencecenter.asterisk.input.AsteriskInputHandler;
import nl.esciencecenter.asterisk.interfaces.SceneStorage;
import nl.esciencecenter.asterisk.interfaces.TimedPlayer;
import nl.esciencecenter.asterisk.interfaces.VisualScene;
import nl.esciencecenter.esight.datastructures.FBO;
import nl.esciencecenter.esight.datastructures.IntPBO;
import nl.esciencecenter.esight.exceptions.UninitializedException;
import nl.esciencecenter.esight.math.Color4;
import nl.esciencecenter.esight.math.MatF4;
import nl.esciencecenter.esight.math.MatrixFMath;
import nl.esciencecenter.esight.math.Point4;
import nl.esciencecenter.esight.math.VecF3;
import nl.esciencecenter.esight.math.VecF4;
import nl.esciencecenter.esight.models.Axis;
import nl.esciencecenter.esight.models.Model;
import nl.esciencecenter.esight.models.Quad;
import nl.esciencecenter.esight.noise.Perlin3D;
import nl.esciencecenter.esight.shaders.ShaderProgram;
import nl.esciencecenter.esight.shaders.ShaderProgramLoader;
import nl.esciencecenter.esight.text.MultiColorText;
import nl.esciencecenter.esight.text.jogampExperimental.Font;
import nl.esciencecenter.esight.text.jogampExperimental.FontFactory;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/* Copyright [2013] [Netherlands eScience Center]
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * @author Maarten van Meersbergen <m.van.meersbergen@esciencecenter.nl>
 * 
 */
public class AsteriskGLEventListener implements GLEventListener {
    private final static Logger logger = LoggerFactory.getLogger(AsteriskGLEventListener.class);

    private ShaderProgram animatedTurbulenceShader, pplShader, starHaloShader, axesShader, pointGasShader,
            postprocessShader, gaussianBlurShader, textShader;

    private FBO starHaloFBO, pointGasFBO, sphereFBO, starFBO, axesFBO, hudFBO;

    private Quad FSQ_postprocess, FSQ_blur;
    private Model xAxis, yAxis, zAxis;

    private final int fontSize = 30;

    private MultiColorText frameNumberText;
    private Perlin3D noiseTex;

    private float offset = 0;

    private final AsteriskSettings settings = AsteriskSettings.getInstance();

    private GlueSceneDescription requestedScene = null;

    private SceneStorage sceneStore;

    private GlueTimedPlayer timer;

    private IntPBO finalPBO;

    private VisualScene oldScene;

    private final AsteriskInputHandler inputHandler;

    protected final ShaderProgramLoader loader;
    protected int canvasWidth, canvasHeight;

    protected int fontSet = FontFactory.UBUNTU;
    protected Font font;
    protected final float radius = 1.0f;
    protected final float ftheta = 0.0f;
    protected final float phi = 0.0f;

    protected final float fovy;
    private float aspect;
    protected final float zNear;
    protected final float zFar;

    private final boolean test = false;

    int lastTime;

    public AsteriskGLEventListener(AsteriskInputHandler inputHandler) {
        this.loader = new ShaderProgramLoader();
        this.inputHandler = inputHandler;
        this.font = FontFactory.get(fontSet).getDefault();

        AsteriskSettings settings = AsteriskSettings.getInstance();

        this.fovy = settings.getFieldOfView();
        this.zNear = settings.getZNear();
        this.zFar = settings.getZFar();

    }

    public static void contextOn(GLAutoDrawable drawable) {
        try {
            final int status = drawable.getContext().makeCurrent();
            if ((status != GLContext.CONTEXT_CURRENT) && (status != GLContext.CONTEXT_CURRENT_NEW)) {
                System.err.println("Error swapping context to onscreen.");
            }
        } catch (final GLException e) {
            System.err.println("Exception while swapping context to onscreen.");
            e.printStackTrace();
        }
    }

    public static void contextOff(GLAutoDrawable drawable) {
        // Release the context.
        try {
            drawable.getContext().release();
        } catch (final GLException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void display(GLAutoDrawable drawable) {
        GlueTimedPlayer timer = AsteriskInterfaceWindow.getTimer();

        if (timer.isInitialized()) {
            this.timer = timer;

            contextOn(drawable);

            final GL3 gl = drawable.getContext().getGL().getGL3();
            gl.glViewport(0, 0, canvasWidth, canvasHeight);

            GlueSceneDescription currentDescription = settings.getCurrentDescription();

            sceneStore = timer.getSceneStorage();
            sceneStore.init(gl);

            if (currentDescription != requestedScene) {
                sceneStore.requestNewConfiguration(currentDescription);
                requestedScene = currentDescription;
            }

            VisualScene newScene = sceneStore.getScene();
            if (newScene != null) {
                ((GlueScene) newScene).init(gl);
                displayContext(newScene);

                if (oldScene != null && oldScene != newScene) {
                    oldScene.dispose(gl);

                    oldScene = newScene;
                }
            } else {
                logger.debug("Scene is null");
            }

            if (timer.isScreenshotNeeded()) {
                finalPBO.makeScreenshotPNG(gl, timer.getScreenshotFileName());

                timer.setScreenshotNeeded(false);
            }

            if (settings.isDoOrbit()) {
                if (settings.isOrbitLinkedToPlayback()) {
                    if (timer.getFrameNumber() != lastTime) {
                        VecF3 rotation = inputHandler.getRotation();
                        VecF3 change = (new VecF3(settings.getXOrbitSpeed(), settings.getYOrbitSpeed(),
                                settings.getZOrbitSpeed())).mul(settings.getOrbitSpeed());
                        rotation = rotation.add(change);
                        inputHandler.setRotation(rotation);
                    }
                    lastTime = timer.getFrameNumber();
                } else {
                    VecF3 rotation = inputHandler.getRotation();
                    VecF3 change = (new VecF3(settings.getXOrbitSpeed(), settings.getYOrbitSpeed(),
                            settings.getZOrbitSpeed())).mul(settings.getOrbitSpeed());
                    rotation = rotation.add(change);
                    inputHandler.setRotation(rotation);
                }
            }

            contextOff(drawable);
        }
    }

    private synchronized void displayContext(VisualScene newScene) {
        final GL3 gl = GLContext.getCurrentGL().getGL3();

        gl.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT);

        final Point4 eye = new Point4((float) (radius * Math.sin(ftheta) * Math.cos(phi)), (float) (radius
                * Math.sin(ftheta) * Math.sin(phi)), (float) (radius * Math.cos(ftheta)), 1.0f);
        final Point4 at = new Point4(0.0f, 0.0f, 0.0f, 1.0f);
        final VecF4 up = new VecF4(0.0f, 1.0f, 0.0f, 0.0f);

        if (settings.getStereo()) {
            MatF4 mv = MatrixFMath.lookAt(eye, at, up);
            mv = mv.mul(MatrixFMath.translate(new VecF3(0f, 0f, inputHandler.getViewDist())));
            MatF4 mv2 = mv.clone();

            if (!settings.getStereoSwitched()) {
                gl.glDrawBuffer(GL2GL3.GL_BACK_LEFT);
            } else {
                gl.glDrawBuffer(GL2GL3.GL_BACK_RIGHT);
            }
            mv = mv.mul(MatrixFMath.translate(new VecF3(-.5f * settings.getStereoOcularDistance(), 0f, 0f)));
            mv = mv.mul(MatrixFMath.rotationX(inputHandler.getRotation().get(0)));
            mv = mv.mul(MatrixFMath.rotationY(inputHandler.getRotation().get(1)));
            mv = mv.mul(MatrixFMath.rotationZ(inputHandler.getRotation().get(2)));
            mv = mv.mul(MatrixFMath.translate(inputHandler.getTranslation()));

            renderScene(gl, mv.clone(), newScene);

            // try {
            // renderHUDText(gl, newScene, mv.clone());
            // } catch (final UninitializedException e) {
            // e.printStackTrace();
            // }
            renderTexturesToScreen(gl);

            if (!settings.getStereoSwitched()) {
                gl.glDrawBuffer(GL2GL3.GL_BACK_RIGHT);
            } else {
                gl.glDrawBuffer(GL2GL3.GL_BACK_LEFT);
            }
            mv2 = mv2.mul(MatrixFMath.translate(new VecF3(.5f * settings.getStereoOcularDistance(), 0f, 0f)));
            mv2 = mv2.mul(MatrixFMath.rotationX(inputHandler.getRotation().get(0)));
            mv2 = mv2.mul(MatrixFMath.rotationY(inputHandler.getRotation().get(1)));
            mv2 = mv2.mul(MatrixFMath.rotationZ(inputHandler.getRotation().get(2)));
            mv2 = mv2.mul(MatrixFMath.translate(inputHandler.getTranslation()));

            renderScene(gl, mv2.clone(), newScene);

            // try {
            // renderHUDText(gl, newScene, mv2.clone());
            // } catch (final UninitializedException e) {
            // e.printStackTrace();
            // }
            renderTexturesToScreen(gl);
        } else {
            MatF4 mv = MatrixFMath.lookAt(eye, at, up);
            mv = mv.mul(MatrixFMath.translate(new VecF3(0f, 0f, inputHandler.getViewDist())));

            mv = mv.mul(MatrixFMath.translate(inputHandler.getTranslation()));
            mv = mv.mul(MatrixFMath.rotationX(inputHandler.getRotation().get(0)));
            mv = mv.mul(MatrixFMath.rotationY(inputHandler.getRotation().get(1)));
            mv = mv.mul(MatrixFMath.rotationZ(inputHandler.getRotation().get(2)));

            renderScene(gl, mv.clone(), newScene);

            try {
                renderHUDText(gl, newScene, mv.clone());
            } catch (final UninitializedException e) {
                e.printStackTrace();
            }
            renderTexturesToScreen(gl);
        }
    }

    private void renderScene(GL3 gl, MatF4 mv, VisualScene newScene) {
        if (settings.getGasInvertedBackgroundColor()) {
            gl.glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        } else {
            gl.glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
        }
        gl.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT);

        try {
            renderPointGas(gl, mv.clone(), newScene);
            renderSpheres(gl, mv, newScene);
            renderStars(gl, mv.clone(), newScene);
            renderStarHalos(gl, mv.clone(), newScene);
            renderAxes(gl, mv.clone());
        } catch (final UninitializedException e) {
            e.printStackTrace();
        }
    }

    private void renderPointGas(GL3 gl, MatF4 mv, VisualScene newScene) throws UninitializedException {
        pointGasFBO.bind(gl);
        gl.glClear(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);

        gl.glDisable(GL.GL_DEPTH_TEST);

        final MatF4 p = MatrixFMath.perspective(fovy, aspect, zNear, zFar);
        pointGasShader.setUniformMatrix("PMatrix", p);
        pointGasShader.setUniform("PointSizeCameraDistDependant",
                settings.isPointgasSizeDependantOnCameraDistance() ? 1 : 0);
        pointGasShader.setUniform("PointSizeMultiplier", (settings.getPointGasPointSizeSetting()));
        pointGasShader.setUniform("CameraDistance", (inputHandler.getViewDist()));

        newScene.drawGasPointCloud(gl, pointGasShader, mv);

        gl.glEnable(GL.GL_DEPTH_TEST);

        pointGasFBO.unBind(gl);

        blur(gl, pointGasFBO, FSQ_blur, settings.getPointGasBlurPassSetting(), settings.getPointGasBlurTypeSetting(),
                settings.getPointGasBlurSizeSetting());
    }

    private void renderSpheres(GL3 gl, MatF4 mv, VisualScene newScene) throws UninitializedException {
        sphereFBO.bind(gl);
        gl.glClear(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);

        final MatF4 p = MatrixFMath.perspective(fovy, aspect, zNear, zFar);
        pplShader.setUniformMatrix("PMatrix", p);
        pplShader.setUniformMatrix("SMatrix", MatrixFMath.scale(settings.getParticleSizeMultiplier() * .1f));
        pplShader.setUniformMatrix("MVMatrix", mv);

        newScene.drawSpheres(gl, pplShader, mv);

        sphereFBO.unBind(gl);
    }

    private void renderStars(GL3 gl, MatF4 mv, VisualScene newScene) throws UninitializedException {
        starFBO.bind(gl);
        gl.glClear(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);

        noiseTex.use(gl);
        animatedTurbulenceShader.setUniform("Noise", noiseTex.getMultitexNumber());

        final MatF4 p = MatrixFMath.perspective(fovy, aspect, zNear, zFar);
        animatedTurbulenceShader.setUniformMatrix("PMatrix", p);
        animatedTurbulenceShader.setUniformMatrix("SMatrix",
                MatrixFMath.scale(settings.getParticleSizeMultiplier() * .1f));
        // animatedTurbulenceShader.setUniformMatrix("MVMatrix", mv);
        animatedTurbulenceShader.setUniform("Offset", offset);

        offset += .001f;

        newScene.drawStars(gl, animatedTurbulenceShader, mv);

        starFBO.unBind(gl);
    }

    private void renderStarHalos(GL3 gl, MatF4 mv, VisualScene newScene) throws UninitializedException {
        starHaloFBO.bind(gl);
        gl.glClear(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);

        final MatF4 p = MatrixFMath.perspective(fovy, aspect, zNear, zFar);
        starHaloShader.setUniformMatrix("PMatrix", p);
        starHaloShader.setUniformMatrix("SMatrix", MatrixFMath.scale((settings.getParticleSizeMultiplier() * .2f)));

        newScene.drawStars(gl, starHaloShader, mv);

        starHaloFBO.unBind(gl);
        blur(gl, starHaloFBO, FSQ_blur, settings.getStarHaloBlurPassSetting(), settings.getStarHaloBlurTypeSetting(),
                settings.getStarHaloBlurSizeSetting());
    }

    private void renderAxes(GL3 gl, MatF4 mv) throws UninitializedException {
        axesFBO.bind(gl);
        gl.glClear(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);

        final MatF4 p = MatrixFMath.perspective(fovy, aspect, zNear, zFar);
        axesShader.setUniformMatrix("PMatrix", p);
        axesShader.setUniformMatrix("MVMatrix", mv);

        axesShader.setUniformVector("Color", new VecF4(1f, 0f, 0f, 1f));
        axesShader.use(gl);
        xAxis.draw(gl, axesShader);

        axesShader.setUniformVector("Color", new VecF4(0f, 1f, 0f, 1f));
        axesShader.use(gl);
        yAxis.draw(gl, axesShader);

        axesShader.setUniformVector("Color", new VecF4(0f, 0f, 1f, 1f));
        axesShader.use(gl);
        zAxis.draw(gl, axesShader);

        axesFBO.unBind(gl);
    }

    private void renderHUDText(GL3 gl, VisualScene newScene, MatF4 mv) throws UninitializedException {
        // SceneDescription currentDesc = newScene.getDescription();
        String frameNumberString = String.format("%03d : ", timer.getFrameNumber());

        frameNumberText.setString(gl, frameNumberString + newScene.getDescriptionString(), Color4.white, fontSize);

        // String min = Float.toString(currentDesc.getLowerBound());
        // String max = Float.toString(currentDesc.getUpperBound());
        // legendTextmin.setString(gl, min, Color4.white, fontSize);
        // legendTextmax.setString(gl, max, Color4.white, fontSize);

        hudFBO.bind(gl);
        gl.glClear(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);

        // Draw Legend
        // ByteBufferTexture legendTexture = new
        // ByteBufferTexture(GL3.GL_TEXTURE9, timer.getSceneStorage()
        // .getLegendImage(), 1, 500);
        // legendTexture.init(gl);
        //
        // legendProgram.setUniform("texture_map",
        // legendTexture.getMultitexNumber());
        // legendProgram.setUniformMatrix("PMatrix", new MatF4());
        // legendProgram.setUniformMatrix("MVMatrix", new MatF4());
        //
        // legendProgram.use(gl);
        // legendModel.draw(gl, legendProgram);
        //
        // legendTexture.delete(gl);
        //
        // // Draw legend text
        // int textLength = legendTextmin.toString().length() * fontSize;
        // legendTextmin.draw(gl, textShader, canvasWidth, canvasHeight, 2 *
        // canvasWidth - textLength - 100,
        // .2f * canvasHeight);
        //
        // textLength = legendTextmax.toString().length() * fontSize;
        // legendTextmax.draw(gl, textShader, canvasWidth, canvasHeight, 2 *
        // canvasWidth - textLength - 100,
        // 1.75f * canvasHeight);

        frameNumberText.draw(gl, textShader, canvasWidth, canvasHeight, 30f, 30f);

        hudFBO.unBind(gl);
    }

    private void renderTexturesToScreen(GL3 gl) {
        gl.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT);

        postprocessShader.setUniform("axesTexture", axesFBO.getTexture().getMultitexNumber());
        postprocessShader.setUniform("pointGasTexture", pointGasFBO.getTexture().getMultitexNumber());
        postprocessShader.setUniform("sphereTexture", sphereFBO.getTexture().getMultitexNumber());
        postprocessShader.setUniform("starTexture", starFBO.getTexture().getMultitexNumber());
        postprocessShader.setUniform("starHaloTexture", starHaloFBO.getTexture().getMultitexNumber());
        postprocessShader.setUniform("hudTexture", hudFBO.getTexture().getMultitexNumber());

        postprocessShader.setUniform("sphereBrightness", settings.getPostprocessingSphereBrightness());
        postprocessShader.setUniform("starBrightness", settings.getPostprocessingStarBrightness());
        postprocessShader.setUniform("starHaloBrightness", settings.getPostprocessingStarHaloBrightness());
        postprocessShader.setUniform("pointGasBrightness", settings.getPostprocessingPointGasBrightness());
        postprocessShader.setUniform("axesBrightness", settings.getPostprocessingAxesBrightness());
        postprocessShader.setUniform("hudBrightness", settings.getPostprocessingHudBrightness());
        postprocessShader.setUniform("overallBrightness", settings.getPostprocessingOverallBrightness());

        postprocessShader.setUniformMatrix("MVMatrix", new MatF4());
        postprocessShader.setUniformMatrix("PMatrix", new MatF4());

        postprocessShader.setUniform("scrWidth", canvasWidth);
        postprocessShader.setUniform("scrHeight", canvasHeight);

        try {
            postprocessShader.use(gl);

            FSQ_blur.draw(gl, postprocessShader);
        } catch (final UninitializedException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void reshape(GLAutoDrawable drawable, int x, int y, int w, int h) {
        contextOn(drawable);

        final GL3 gl = GLContext.getCurrentGL().getGL3();

        canvasWidth = GLContext.getCurrent().getGLDrawable().getWidth();
        canvasHeight = GLContext.getCurrent().getGLDrawable().getHeight();

        aspect = (float) canvasWidth / (float) canvasHeight;

        starFBO.delete(gl);
        starHaloFBO.delete(gl);
        pointGasFBO.delete(gl);
        axesFBO.delete(gl);
        hudFBO.delete(gl);
        sphereFBO.delete(gl);

        starFBO = new FBO(canvasWidth, canvasHeight, GL.GL_TEXTURE1);
        starHaloFBO = new FBO(canvasWidth, canvasHeight, GL.GL_TEXTURE2);
        pointGasFBO = new FBO(canvasWidth, canvasHeight, GL.GL_TEXTURE3);
        axesFBO = new FBO(canvasWidth, canvasHeight, GL.GL_TEXTURE4);
        hudFBO = new FBO(canvasWidth, canvasHeight, GL.GL_TEXTURE5);
        sphereFBO = new FBO(canvasWidth, canvasHeight, GL.GL_TEXTURE6);

        starFBO.init(gl);
        starHaloFBO.init(gl);
        pointGasFBO.init(gl);
        axesFBO.init(gl);
        hudFBO.init(gl);
        sphereFBO.init(gl);

        finalPBO.delete(gl);
        finalPBO = new IntPBO(canvasWidth, canvasHeight);
        finalPBO.init(gl);

        contextOff(drawable);
    }

    private void blur(GL3 gl, FBO target, Quad fullScreenQuad, int passes, int blurType, float blurSize) {
        gaussianBlurShader.setUniform("Texture", target.getTexture().getMultitexNumber());

        gaussianBlurShader.setUniformMatrix("PMatrix", new MatF4());
        gaussianBlurShader.setUniformMatrix("MVMatrix", new MatF4());

        gaussianBlurShader.setUniform("blurType", blurType);
        gaussianBlurShader.setUniform("blurSize", blurSize);

        gaussianBlurShader.setUniform("scrWidth", target.getTexture().getWidth());
        gaussianBlurShader.setUniform("scrHeight", target.getTexture().getHeight());

        gaussianBlurShader.setUniform("Alpha", 1f);
        gaussianBlurShader.setUniform("ColorMultiplier", 1.25f);

        // gaussianBlurShader.setUniform("blurDirection", 0);

        gaussianBlurShader.setUniform("NumPixelsPerSide", 2f);
        gaussianBlurShader.setUniform("Sigma", 2f);

        try {
            target.bind(gl);
            for (int i = 0; i < passes; i++) {
                gaussianBlurShader.setUniform("blurDirection", 0);
                gaussianBlurShader.use(gl);
                fullScreenQuad.draw(gl, gaussianBlurShader);

                gaussianBlurShader.setUniform("blurDirection", 1);
                gaussianBlurShader.use(gl);
                fullScreenQuad.draw(gl, gaussianBlurShader);
            }

            target.unBind(gl);
        } catch (final UninitializedException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void dispose(GLAutoDrawable drawable) {
        contextOn(drawable);

        final GL3 gl = drawable.getGL().getGL3();

        noiseTex.delete(gl);

        starFBO.delete(gl);
        starHaloFBO.delete(gl);
        axesFBO.delete(gl);
        hudFBO.delete(gl);
        sphereFBO.delete(gl);

        finalPBO.delete(gl);

        try {
            loader.cleanup(gl);
        } catch (UninitializedException e) {
            e.printStackTrace();
        }

        contextOff(drawable);
    }

    @Override
    public void init(GLAutoDrawable drawable) {
        contextOn(drawable);

        final GL3 gl = GLContext.getCurrentGL().getGL3();

        canvasWidth = GLContext.getCurrent().getGLDrawable().getWidth();
        canvasHeight = GLContext.getCurrent().getGLDrawable().getHeight();

        aspect = (float) canvasWidth / (float) canvasHeight;

        // Anti-Aliasing
        gl.glEnable(GL3.GL_LINE_SMOOTH);
        gl.glHint(GL3.GL_LINE_SMOOTH_HINT, GL3.GL_NICEST);
        gl.glEnable(GL3.GL_POLYGON_SMOOTH);
        gl.glHint(GL3.GL_POLYGON_SMOOTH_HINT, GL3.GL_NICEST);

        // Depth testing
        gl.glEnable(GL3.GL_DEPTH_TEST);
        gl.glDepthFunc(GL3.GL_LEQUAL);
        gl.glClearDepth(1.0f);

        // Culling
        gl.glEnable(GL3.GL_CULL_FACE);
        gl.glCullFace(GL3.GL_BACK);

        // Enable Blending (needed for both Transparency and
        // Anti-Aliasing
        gl.glBlendFunc(GL3.GL_SRC_ALPHA, GL3.GL_ONE_MINUS_SRC_ALPHA);
        gl.glEnable(GL3.GL_BLEND);

        // Enable Vertical Sync
        gl.setSwapInterval(1);

        // Set black background
        gl.glClearColor(0f, 0f, 0f, 0f);

        // enable setting of point sizes
        gl.glEnable(GL3.GL_PROGRAM_POINT_SIZE);

        // Load and compile shaders, then use program.
        try {
            animatedTurbulenceShader = loader.createProgram(gl, "animatedTurbulence", new File(
                    "shaders/vs_sunsurface.vp"), new File("shaders/fs_animatedTurbulence.fp"));
            pplShader = loader.createProgram(gl, "ppl", new File("shaders/vs_ppl.vp"), new File("shaders/fs_ppl.fp"));
            starHaloShader = loader.createProgram(gl, "starHalo", new File("shaders/vs_starHalo.vp"), new File(
                    "shaders/fs_starHalo.fp"));
            axesShader = loader.createProgram(gl, "axes", new File("shaders/vs_axes.vp"),
                    new File("shaders/fs_axes.fp"));
            pointGasShader = loader.createProgram(gl, "gas", new File("shaders/vs_gas.vp"), new File(
                    "shaders/fs_gas.fp"));
            textShader = loader.createProgram(gl, "text", new File("shaders/vs_multiColorTextShader.vp"), new File(
                    "shaders/fs_multiColorTextShader.fp"));
            postprocessShader = loader.createProgram(gl, "postprocess", new File("shaders/vs_postprocess.vp"),
                    new File("shaders/fs_postprocess.fp"));
            gaussianBlurShader = loader.createProgram(gl, "gaussianBlur", new File("shaders/vs_postprocess.vp"),
                    new File("shaders/fs_gaussian_blur.fp"));
        } catch (final Exception e) {
            System.err.println(e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }

        // AXES
        xAxis = new Axis(new VecF3(-1f, 0f, 0f), new VecF3(1f, 0f, 0f), .1f, .02f);
        xAxis.init(gl);
        yAxis = new Axis(new VecF3(0f, -1f, 0f), new VecF3(0f, 1f, 0f), .1f, .02f);
        yAxis.init(gl);
        zAxis = new Axis(new VecF3(0f, 0f, -1f), new VecF3(0f, 0f, 1f), .1f, .02f);
        zAxis.init(gl);

        // TEXT
        String text = "Frame: " + AsteriskInterfaceWindow.getTimer().getFrameNumber();
        frameNumberText = new MultiColorText(gl, font, text, Color4.white, fontSize);
        // frameNumberText.setString(gl, text, Color4.white, fontSize);

        frameNumberText.init(gl);

        // FULL SCREEN QUADS
        FSQ_postprocess = new Quad(2, 2, new VecF3(0, 0, 0.1f));
        FSQ_postprocess.init(gl);

        FSQ_blur = new Quad(2, 2, new VecF3(0, 0, 0.1f));
        FSQ_blur.init(gl);

        // TEXTURES
        noiseTex = new Perlin3D(GL.GL_TEXTURE0, 128, 128, 128);
        noiseTex.init(gl);

        // Full screen textures (for post processing) done with FBO's
        starFBO = new FBO(canvasWidth, canvasHeight, GL.GL_TEXTURE1);
        starHaloFBO = new FBO(canvasWidth, canvasHeight, GL.GL_TEXTURE2);
        pointGasFBO = new FBO(canvasWidth, canvasHeight, GL.GL_TEXTURE3);
        axesFBO = new FBO(canvasWidth, canvasHeight, GL.GL_TEXTURE4);
        hudFBO = new FBO(canvasWidth, canvasHeight, GL.GL_TEXTURE5);
        sphereFBO = new FBO(canvasWidth, canvasHeight, GL.GL_TEXTURE6);

        starFBO.init(gl);
        starHaloFBO.init(gl);
        pointGasFBO.init(gl);
        axesFBO.init(gl);
        hudFBO.init(gl);
        sphereFBO.init(gl);

        finalPBO = new IntPBO(canvasWidth, canvasHeight);
        finalPBO.init(gl);

        if (test) {
            int POINTS = 10000000;

            PointGas[] pGas1 = new PointGas[POINTS];
            int tenPercent = (int) (0.1f * pGas1.length);
            long startTime = System.currentTimeMillis();

            for (int i = 0; i < pGas1.length; i++) {
                float[] coordinates = new float[] { randBound(-1f, 2f), randBound(-1f, 2f), randBound(-1f, 2f) };

                float[] color = new float[] {
                        // 1f, 1f, 1f, 1f
                        (float) Math.random(), (float) Math.random(), (float) Math.random(), (float) Math.random() };

                pGas1[i] = new PointGas(i, coordinates, color);

                if (i % tenPercent == 0) {
                    long stopTime = System.currentTimeMillis();
                    long timePassed = stopTime - startTime;

                    System.out.println("Point cloud generation progress: "
                            + (int) (((float) i / (float) pGas1.length) * 100) + "% ... Time passed: "
                            + TimeUnit.MILLISECONDS.toSeconds(timePassed) + " seconds.");

                }
            }

            Snapshot scene1 = new Snapshot("willekeurig", null, null, null, pGas1);

            TimedPlayer timer = AsteriskInterfaceWindow.getTimer();
            int sceneNumber = ((GlueTimedPlayer) timer).addScene(scene1);

            if (!timer.isInitialized()) {
                ((GlueTimedPlayer) timer).init();

                new Thread(timer).start();
            }
        }

        contextOff(drawable);
    }

    private static float randBound(float lower, float delta) {
        return ((float) Math.random() * delta + lower);
    }

    public AsteriskInputHandler getInputHandler() {
        return inputHandler;
    }

    public TimedPlayer getTimer() {
        return timer;
    }
}
