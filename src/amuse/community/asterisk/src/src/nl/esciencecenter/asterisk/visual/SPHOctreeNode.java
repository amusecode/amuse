package nl.esciencecenter.asterisk.visual;

import java.util.ArrayList;

import javax.media.opengl.GL3;

import nl.esciencecenter.asterisk.AsteriskSettings;
import nl.esciencecenter.asterisk.data.GlueSceneDescription;
import nl.esciencecenter.esight.exceptions.UninitializedException;
import nl.esciencecenter.esight.input.InputHandler;
import nl.esciencecenter.esight.math.MatF4;
import nl.esciencecenter.esight.math.MatrixFMath;
import nl.esciencecenter.esight.math.VecF3;
import nl.esciencecenter.esight.math.VecF4;
import nl.esciencecenter.esight.models.Model;
import nl.esciencecenter.esight.shaders.ShaderProgram;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class SPHOctreeNode {
    private final static AsteriskSettings settings = AsteriskSettings
            .getInstance();
    private final static Logger logger = LoggerFactory
            .getLogger(SPHOctreeNode.class);

    public static enum Octant {
        PPP, PPN, PNP, PNN, NPP, NPN, NNP, NNN
    };

    private class SPHOctreeElement {
        private final VecF3 center;
        private final VecF4 color;

        public SPHOctreeElement(VecF3 center, VecF4 color) {
            this.center = center;
            this.color = color;
        }

        public VecF3 getCenter() {
            return center;
        }

        public VecF4 getColor() {
            return color;
        }
    }

    private final ArrayList<SPHOctreeElement> elements;
    private final VecF3 center;
    private final float cubeSize;
    private final int depth;
    private final Model baseModel;
    private final MatF4 TMatrix;
    private final float modelScale;

    private SPHOctreeNode ppp, ppn, pnp, npp, pnn, npn, nnp, nnn;
    private int childCounter;
    private boolean subdivided = false;
    private boolean initialized = false;
    private VecF4 color;

    private final GlueSceneDescription description;

    public SPHOctreeNode(Model baseModel, int depth, VecF3 center, float size,
            GlueSceneDescription description) {
        this.baseModel = baseModel;
        this.depth = depth;
        this.center = center;
        this.cubeSize = size;
        this.modelScale = size * 3f;
        this.elements = new ArrayList<SPHOctreeElement>();
        this.childCounter = 0;
        this.description = description;

        this.TMatrix = MatrixFMath.translate(center);
    }

    public void addElement(VecF3 coordinates, VecF4 color) {
        SPHOctreeElement element = new SPHOctreeElement(coordinates, color);
        addElement(element);
    }

    public void addElement(SPHOctreeElement element) {
        final VecF3 location = element.getCenter();

        if ((location.get(0) > (center.get(0) - cubeSize))
                && (location.get(1) > (center.get(1) - cubeSize))
                && (location.get(2) > (center.get(2) - cubeSize))
                && (location.get(0) < (center.get(0) + cubeSize))
                && (location.get(1) < (center.get(1) + cubeSize))
                && (location.get(2) < (center.get(2) + cubeSize))) {

            if (subdivided) {
                addElementSubdivided(element);
            } else {
                if (childCounter > AsteriskSettings
                        .getMaxOctreeElementsPerNode()) {
                    if (depth < settings.getMaxOctreeDepth()) {
                        subDiv();
                        childCounter = 0;
                    } else {
                        logger.warn("Octree is maximally divided, please change MAX_OCTREE_DEPTH constant if you are certain.");
                    }
                } else {
                    elements.add(element);

                    childCounter++;
                }
            }
        }
    }

    public void addElementSubdivided(SPHOctreeElement element) {
        VecF3 location = element.getCenter();
        if (location.get(0) < center.get(0)) {
            if (location.get(1) < center.get(1)) {
                if (location.get(2) < center.get(2)) {
                    nnn.addElement(element);
                } else {
                    nnp.addElement(element);
                }
            } else {
                if (location.get(2) < center.get(2)) {
                    npn.addElement(element);
                } else {
                    npp.addElement(element);
                }
            }
        } else {
            if (location.get(1) < center.get(1)) {
                if (location.get(2) < center.get(2)) {
                    pnn.addElement(element);
                } else {
                    pnp.addElement(element);
                }
            } else {
                if (location.get(2) < center.get(2)) {
                    ppn.addElement(element);
                } else {
                    ppp.addElement(element);
                }
            }
        }
    }

    protected void subDiv() {
        float childSize = cubeSize / 2f;

        VecF3 pppCenter = center
                .add(new VecF3(childSize, childSize, childSize));
        VecF3 ppnCenter = center
                .add(new VecF3(childSize, childSize, -childSize));
        VecF3 pnpCenter = center
                .add(new VecF3(childSize, -childSize, childSize));
        VecF3 nppCenter = center
                .add(new VecF3(-childSize, childSize, childSize));
        VecF3 pnnCenter = center.add(new VecF3(childSize, -childSize,
                -childSize));
        VecF3 npnCenter = center.add(new VecF3(-childSize, childSize,
                -childSize));
        VecF3 nnpCenter = center.add(new VecF3(-childSize, -childSize,
                childSize));
        VecF3 nnnCenter = center.add(new VecF3(-childSize, -childSize,
                -childSize));

        ppp = new SPHOctreeNode(baseModel, depth + 1, pppCenter, childSize,
                description);
        ppn = new SPHOctreeNode(baseModel, depth + 1, ppnCenter, childSize,
                description);
        pnp = new SPHOctreeNode(baseModel, depth + 1, pnpCenter, childSize,
                description);
        npp = new SPHOctreeNode(baseModel, depth + 1, nppCenter, childSize,
                description);
        pnn = new SPHOctreeNode(baseModel, depth + 1, pnnCenter, childSize,
                description);
        npn = new SPHOctreeNode(baseModel, depth + 1, npnCenter, childSize,
                description);
        nnp = new SPHOctreeNode(baseModel, depth + 1, nnpCenter, childSize,
                description);
        nnn = new SPHOctreeNode(baseModel, depth + 1, nnnCenter, childSize,
                description);

        for (SPHOctreeElement element : elements) {
            addElementSubdivided(element);
        }

        elements.clear();

        subdivided = true;
    }

    public void draw(GL3 gl, ShaderProgram program) {
        if (subdivided) {
            draw_sorted(gl, program);
        } else {

            // System.err.println("alpha " + density);
            if (color.get(3) > 0.01f) {
                program.setUniform("node_scale", modelScale);
                program.setUniformVector("Color", color);
                program.setUniformMatrix("node_MVmatrix", TMatrix);

                program.setUniformMatrix("node_NormalMatrix",
                        MatrixFMath.getNormalMatrix(TMatrix));

                try {
                    program.use(gl);
                } catch (UninitializedException e) {
                    e.printStackTrace();
                }

                baseModel.draw(gl, program);
            }
        }
    }

    protected void draw_unsorted(GL3 gl, ShaderProgram program) {
        nnn.draw(gl, program);
        pnn.draw(gl, program);
        npn.draw(gl, program);
        nnp.draw(gl, program);
        ppn.draw(gl, program);
        npp.draw(gl, program);
        pnp.draw(gl, program);
        ppp.draw(gl, program);
    }

    protected void draw_sorted(GL3 gl, ShaderProgram program) {
        InputHandler inputHandler = InputHandler.getInstance();

        if (inputHandler.getCurrentViewOctant() == InputHandler.octants.NNN) {
            ppp.draw(gl, program);

            npp.draw(gl, program);
            pnp.draw(gl, program);
            ppn.draw(gl, program);

            nnp.draw(gl, program);
            pnn.draw(gl, program);
            npn.draw(gl, program);

            nnn.draw(gl, program);
        } else if (inputHandler.getCurrentViewOctant() == InputHandler.octants.NNP) {
            ppn.draw(gl, program);

            npn.draw(gl, program);
            pnn.draw(gl, program);
            ppp.draw(gl, program);

            nnn.draw(gl, program);
            pnp.draw(gl, program);
            npp.draw(gl, program);

            nnp.draw(gl, program);
        } else if (inputHandler.getCurrentViewOctant() == InputHandler.octants.NPN) {
            pnp.draw(gl, program);

            nnp.draw(gl, program);
            ppp.draw(gl, program);
            pnn.draw(gl, program);

            npp.draw(gl, program);
            ppn.draw(gl, program);
            nnn.draw(gl, program);

            npn.draw(gl, program);
        } else if (inputHandler.getCurrentViewOctant() == InputHandler.octants.NPP) {
            pnn.draw(gl, program);

            nnn.draw(gl, program);
            ppn.draw(gl, program);
            pnp.draw(gl, program);

            npn.draw(gl, program);
            ppp.draw(gl, program);
            nnp.draw(gl, program);

            npp.draw(gl, program);
        } else if (inputHandler.getCurrentViewOctant() == InputHandler.octants.PNN) {
            npp.draw(gl, program);

            ppp.draw(gl, program);
            nnp.draw(gl, program);
            npn.draw(gl, program);

            pnp.draw(gl, program);
            nnn.draw(gl, program);
            ppn.draw(gl, program);

            pnn.draw(gl, program);
        } else if (inputHandler.getCurrentViewOctant() == InputHandler.octants.PNP) {
            npn.draw(gl, program);

            ppn.draw(gl, program);
            nnn.draw(gl, program);
            npp.draw(gl, program);

            pnn.draw(gl, program);
            nnp.draw(gl, program);
            ppp.draw(gl, program);

            pnp.draw(gl, program);
        } else if (inputHandler.getCurrentViewOctant() == InputHandler.octants.PPN) {
            nnp.draw(gl, program);

            pnp.draw(gl, program);
            npp.draw(gl, program);
            nnn.draw(gl, program);

            ppp.draw(gl, program);
            npn.draw(gl, program);
            pnn.draw(gl, program);

            ppn.draw(gl, program);
        } else if (inputHandler.getCurrentViewOctant() == InputHandler.octants.PPP) {
            nnn.draw(gl, program);

            pnn.draw(gl, program);
            npn.draw(gl, program);
            nnp.draw(gl, program);

            ppn.draw(gl, program);
            npp.draw(gl, program);
            pnp.draw(gl, program);

            ppp.draw(gl, program);
        }
    }

    public void init() {
        if (!initialized) {
            if (subdivided) {
                ppp.init();
                ppn.init();
                pnp.init();
                npp.init();
                pnn.init();
                npn.init();
                nnp.init();
                nnn.init();
            } else {
                VecF4 tmpColor = new VecF4();
                if (elements.size() != 0) {
                    for (SPHOctreeElement element : elements) {
                        tmpColor = tmpColor.add(element.getColor());
                    }
                    color = tmpColor.div(elements.size());

                    float alpha = (elements.size() / AsteriskSettings
                            .getMaxOctreeElementsPerNode())
                            * (depth * depth * depth)
                            / (float) (AsteriskSettings.MAX_OCTREE_DEPTH
                                    * AsteriskSettings.MAX_OCTREE_DEPTH * AsteriskSettings.MAX_OCTREE_DEPTH);
                    color.set(3, alpha * alpha);

                    // System.err.println(color);
                } else {
                    color = tmpColor;
                }
            }

            elements.clear();
            initialized = true;
        }
    }
    // public FloatBuffer gatherTriangles() {
    // if (subdivided) {
    // FloatBuffer pppBuffer = ppp.gatherTriangles();
    // FloatBuffer ppnBuffer = ppn.gatherTriangles();
    // FloatBuffer pnpBuffer = pnp.gatherTriangles();
    // FloatBuffer nppBuffer = npp.gatherTriangles();
    // FloatBuffer pnnBuffer = pnn.gatherTriangles();
    // FloatBuffer npnBuffer = npn.gatherTriangles();
    // FloatBuffer nnpBuffer = nnp.gatherTriangles();
    // FloatBuffer nnnBuffer = nnn.gatherTriangles();
    //
    // int size = pppBuffer.capacity() +
    // ppnBuffer.capacity() +
    // pnpBuffer.capacity() +
    // nppBuffer.capacity() +
    // pnnBuffer.capacity() +
    // npnBuffer.capacity() +
    // nnpBuffer.capacity() +
    // nnnBuffer.capacity();
    //
    // FloatBuffer result = Buffers.newDirectFloatBuffer(size);
    //
    // result.put(pppBuffer);
    // result.put(ppnBuffer);
    // result.put(pnpBuffer);
    // result.put(nppBuffer);
    // result.put(pnnBuffer);
    // result.put(npnBuffer);
    // result.put(nnpBuffer);
    // result.put(nnnBuffer);
    //
    // result.rewind();
    // } else {
    // baseModel
    // }
    // }
}
