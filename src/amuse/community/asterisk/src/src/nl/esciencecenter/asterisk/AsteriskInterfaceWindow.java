package nl.esciencecenter.asterisk;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.FocusTraversalPolicy;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JMenuBar;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.JTabbedPane;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.plaf.basic.BasicSliderUI;

import nl.esciencecenter.asterisk.data.GlueTimedPlayer;
import nl.esciencecenter.asterisk.input.AsteriskInputHandler;
import nl.esciencecenter.esight.math.VecF3;
import nl.esciencecenter.esight.swing.CustomJSlider;
import nl.esciencecenter.esight.swing.GoggleSwing;

public class AsteriskInterfaceWindow extends JPanel {
    private final AsteriskSettings settings = AsteriskSettings.getInstance();

    private static final long serialVersionUID = 1L;

    protected CustomJSlider timeBar;

    protected JFormattedTextField frameCounter;

    public static GlueTimedPlayer timer;

    private final JTabbedPane configPanel;

    private final JPanel starConfig, pointCloudConfig, miscConfig, viewConfig, recordingConfig;

    private final AsteriskInputHandler inputHandler = AsteriskInputHandler.getInstance();

    public class KeyFrame {
        private Component uiElement;
        private final int frameNumber;
        private VecF3 rotation;
        private float viewDist;

        public KeyFrame(int frameNumber) {
            this.frameNumber = frameNumber;
        }

        public Component getUiElement() {
            return uiElement;
        }

        public void setUiElement(Component uiElement) {
            this.uiElement = uiElement;
        }

        public int getFrameNumber() {
            return frameNumber;
        }

        public VecF3 getRotation() {
            return rotation;
        }

        public void setRotation(VecF3 rotation) {
            this.rotation = rotation;
        }

        public float getViewDist() {
            return viewDist;
        }

        public void setViewDist(float viewDist) {
            this.viewDist = viewDist;
        }

        @Override
        public int hashCode() {
            int hash = 1;
            hash = hash * 17 + frameNumber;
            return hash;
        }

        @Override
        public boolean equals(Object other) {
            if (other instanceof KeyFrame && ((KeyFrame) other).hashCode() == this.hashCode()) {
                return true;
            } else {
                return false;
            }
        }
    }

    private final ArrayList<KeyFrame> keyFrames;

    public AsteriskInterfaceWindow() {
        keyFrames = new ArrayList<KeyFrame>();

        setLayout(new BorderLayout(0, 0));

        timeBar = new CustomJSlider(new BasicSliderUI(timeBar));
        timeBar.setValue(0);
        timeBar.setMajorTickSpacing(5);
        timeBar.setMinorTickSpacing(1);
        timeBar.setMaximum(0);
        timeBar.setMinimum(0);
        timeBar.setPaintTicks(true);
        timeBar.setSnapToTicks(true);

        timer = new GlueTimedPlayer(timeBar, frameCounter);
        timer.setScreenshotDirectory(System.getProperty("user.dir") + "/screenshots/");
        settings.setScreenshotPath(System.getProperty("user.dir") + "/screenshots/");

        // Make the menu bar
        final JMenuBar menuBar = new JMenuBar();
        menuBar.setLayout(new BoxLayout(menuBar, BoxLayout.X_AXIS));

        final JMenuBar menuBar2 = new JMenuBar();

        ImageIcon nlescIcon = GoggleSwing.createResizedImageIcon("images/ESCIENCE_logo.jpg", "eScienceCenter Logo",
                200, 20);
        JLabel nlesclogo = new JLabel(nlescIcon);

        menuBar2.add(Box.createHorizontalGlue());
        menuBar2.add(nlesclogo);
        menuBar2.add(Box.createHorizontalGlue());

        Container menuContainer = new Container();
        menuContainer.setLayout(new BoxLayout(menuContainer, BoxLayout.Y_AXIS));

        menuContainer.add(menuBar);
        menuContainer.add(menuBar2);

        add(menuContainer, BorderLayout.NORTH);

        // Make the "media player" panel
        final JPanel bottomPanel = createBottomPanel();

        // Add the tweaks panels
        configPanel = new JTabbedPane();
        add(configPanel, BorderLayout.WEST);
        // configPanel.setLayout(new BoxLayout(configPanel, BoxLayout.Y_AXIS));
        configPanel.setPreferredSize(new Dimension(240, 10));

        miscConfig = new JPanel();
        miscConfig.setLayout(new BoxLayout(miscConfig, BoxLayout.Y_AXIS));
        miscConfig.setMinimumSize(configPanel.getPreferredSize());
        createMiscPanel(miscConfig);
        configPanel.addTab("Misc", miscConfig);

        starConfig = new JPanel();
        starConfig.setLayout(new BoxLayout(starConfig, BoxLayout.Y_AXIS));
        starConfig.setMinimumSize(configPanel.getPreferredSize());
        createStarPanel(starConfig);
        configPanel.addTab("Stars", starConfig);

        pointCloudConfig = new JPanel();
        pointCloudConfig.setLayout(new BoxLayout(pointCloudConfig, BoxLayout.Y_AXIS));
        pointCloudConfig.setMinimumSize(configPanel.getPreferredSize());
        createPointCloudPanel(pointCloudConfig);
        configPanel.addTab("Point Gas", pointCloudConfig);

        viewConfig = new JPanel();
        viewConfig.setLayout(new BoxLayout(viewConfig, BoxLayout.Y_AXIS));
        viewConfig.setMinimumSize(configPanel.getPreferredSize());
        createViewPanel(viewConfig);
        configPanel.addTab("View", viewConfig);

        recordingConfig = new JPanel();
        recordingConfig.setLayout(new BoxLayout(recordingConfig, BoxLayout.Y_AXIS));
        recordingConfig.setMinimumSize(configPanel.getPreferredSize());
        createRecordingPanel(recordingConfig);
        configPanel.addTab("Recording", recordingConfig);

        configPanel.setVisible(true);

        add(bottomPanel, BorderLayout.SOUTH);

        // Read command line file information
        makeTimer();

        setVisible(true);
    }

    private void createRecordingPanel(JPanel targetPanel) {
        ActionListener addKeyFrameListener = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                int frameNumber = settings.getCurrentDescription().getFrameNumber();
                final KeyFrame newKeyFrame = new KeyFrame(frameNumber);

                if (keyFrames.contains(newKeyFrame)) {
                    keyFrames.remove(newKeyFrame);
                }

                ActionListener removeKeyFrameListener = new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        keyFrames.remove(newKeyFrame);
                        recordingConfig.removeAll();
                        createRecordingPanel(recordingConfig);
                        validate();
                        repaint();
                    }
                };

                ArrayList<Component> axesVboxList = new ArrayList<Component>();
                ArrayList<Component> axesHboxList2 = new ArrayList<Component>();

                VecF3 rotation = inputHandler.getRotation().clone();
                newKeyFrame.setRotation(rotation);
                float viewDist = inputHandler.getViewDist();
                newKeyFrame.setViewDist(viewDist);

                axesHboxList2.add(new JLabel("#: " + frameNumber + " Axes: " + (int) rotation.get(0) + " "
                        + (int) rotation.get(1) + " " + (int) rotation.get(2)));
                axesHboxList2.add(Box.createHorizontalGlue());
                JButton removeButton = new JButton(new ImageIcon("images/RemoveIcon15.png"));
                removeButton.addActionListener(removeKeyFrameListener);
                axesHboxList2.add(removeButton);

                axesVboxList.add(GoggleSwing.hBoxedComponents(axesHboxList2, false));
                newKeyFrame.setUiElement(GoggleSwing.vBoxedComponents(axesVboxList, true));

                keyFrames.add(newKeyFrame);

                recordingConfig.removeAll();
                createRecordingPanel(recordingConfig);
                validate();
                repaint();
            }
        };

        ActionListener playSequenceListener = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                timer.startSequence(keyFrames, false);
            }
        };

        ActionListener recordSequenceListener = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                timer.startSequence(keyFrames, true);
            }
        };

        ActionListener clearListListener = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                keyFrames.clear();
                recordingConfig.removeAll();
                createRecordingPanel(recordingConfig);
                validate();
                repaint();
            }
        };

        for (KeyFrame keyFrame : keyFrames) {
            targetPanel.add(keyFrame.getUiElement());
        }

        targetPanel.add(GoggleSwing.buttonBox("", new GoggleSwing.ButtonBoxItem("Add current", addKeyFrameListener),
                new GoggleSwing.ButtonBoxItem("Play Sequence", playSequenceListener), new GoggleSwing.ButtonBoxItem(
                        "Record Sequence", recordSequenceListener), new GoggleSwing.ButtonBoxItem("Clear All",
                        clearListListener)));

    }

    private void createViewPanel(JPanel targetPanel) {
        final ItemListener orbitCheckboxListener = new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent e) {
                final JCheckBox source = (JCheckBox) e.getSource();
                settings.setDoOrbit(source.isSelected());
            }
        };

        targetPanel.add(GoggleSwing.checkboxBox("", new GoggleSwing.CheckBoxItem("Orbit?", settings.isDoOrbit(),
                orbitCheckboxListener)));

        final ItemListener orbitLinkedToPlaybackListener = new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent e) {
                final JCheckBox source = (JCheckBox) e.getSource();
                settings.setOrbitLinkedToPlayback(source.isSelected());
            }
        };

        targetPanel.add(GoggleSwing.checkboxBox("",
                new GoggleSwing.CheckBoxItem("Orbit Linked to Playback?", settings.isOrbitLinkedToPlayback(),
                        orbitLinkedToPlaybackListener)));

        final ChangeListener orbitSpeedSliderListener = new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                final JSlider source = (JSlider) e.getSource();
                if (source.hasFocus()) {
                    settings.setOrbitSpeed(source.getValue() / 10f);
                }
            }
        };
        targetPanel.add(GoggleSwing.sliderBox("Orbit Speed", orbitSpeedSliderListener, 0f, 10f, 1f,
                settings.getOrbitSpeed() * 10f, new JLabel("")));

        JFormattedTextField xOrbit = new JFormattedTextField();
        xOrbit.setValue(new Float(0f));
        xOrbit.setColumns(4);
        xOrbit.setMaximumSize(new Dimension(40, 20));
        xOrbit.addPropertyChangeListener(new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent e) {
                final JFormattedTextField source = (JFormattedTextField) e.getSource();
                settings.setXOrbitSpeed((Float) (source.getValue()));
            }
        });

        JFormattedTextField yOrbit = new JFormattedTextField();
        yOrbit.setValue(new Float(1f));
        yOrbit.setColumns(4);
        yOrbit.setMaximumSize(new Dimension(40, 20));
        yOrbit.addPropertyChangeListener(new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent e) {
                final JFormattedTextField source = (JFormattedTextField) e.getSource();
                settings.setYOrbitSpeed((Float) source.getValue());
            }
        });

        JFormattedTextField zOrbit = new JFormattedTextField();
        zOrbit.setValue(new Float(0f));
        zOrbit.setColumns(4);
        zOrbit.setMaximumSize(new Dimension(40, 20));
        zOrbit.addPropertyChangeListener(new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent e) {
                final JFormattedTextField source = (JFormattedTextField) e.getSource();
                settings.setZOrbitSpeed((Float) source.getValue());
            }
        });

        ArrayList<Component> axesVboxList = new ArrayList<Component>();
        ArrayList<Component> axesHboxList1 = new ArrayList<Component>();
        ArrayList<Component> axesHboxList2 = new ArrayList<Component>();
        axesHboxList1.add(Box.createHorizontalGlue());

        axesHboxList1.add(new JLabel("X:"));
        axesHboxList1.add(xOrbit);
        axesHboxList1.add(GoggleSwing.horizontalStrut(3));
        axesHboxList1.add(new JLabel("Y:"));
        axesHboxList1.add(yOrbit);
        axesHboxList1.add(GoggleSwing.horizontalStrut(3));
        axesHboxList1.add(new JLabel("Z:"));
        axesHboxList1.add(zOrbit);

        axesHboxList1.add(Box.createHorizontalGlue());

        axesHboxList2.add(new JLabel("Orbit direction:"));
        axesHboxList2.add(Box.createHorizontalGlue());

        axesVboxList.add(GoggleSwing.hBoxedComponents(axesHboxList2, false));
        axesVboxList.add(GoggleSwing.hBoxedComponents(axesHboxList1, false));
        targetPanel.add(GoggleSwing.vBoxedComponents(axesVboxList, true));

    }

    private void createMiscPanel(JPanel targetPanel) {
        final ChangeListener overallBrightnessSliderListener = new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                final JSlider source = (JSlider) e.getSource();
                if (source.hasFocus()) {
                    settings.setPostprocessingOverallBrightness(source.getValue());
                }
            }
        };
        targetPanel.add(GoggleSwing.sliderBox("Overall Brightness", overallBrightnessSliderListener,
                (int) (settings.getPostprocessingOverallBrightnessMin()),
                (int) (settings.getPostprocessingOverallBrightnessMax()), 1,
                (int) (settings.getPostprocessingOverallBrightness()), new JLabel("")));

        targetPanel.add(GoggleSwing.verticalStrut(1));

        final ChangeListener axesBrightnessSliderListener = new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                final JSlider source = (JSlider) e.getSource();
                if (source.hasFocus()) {
                    settings.setPostprocessingAxesBrightness(source.getValue());
                }
            }
        };
        targetPanel.add(GoggleSwing.sliderBox("Axes Brightness", axesBrightnessSliderListener,
                (int) (settings.getPostprocessingAxesBrightnessMin()),
                (int) (settings.getPostprocessingAxesBrightnessMax()), 1,
                (int) (settings.getPostprocessingAxesBrightness()), new JLabel("")));

        targetPanel.add(GoggleSwing.verticalStrut(1));

        final ChangeListener hudBrightnessSliderListener = new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                final JSlider source = (JSlider) e.getSource();
                if (source.hasFocus()) {
                    settings.setPostprocessingHudBrightness(source.getValue());
                }
            }
        };
        targetPanel.add(GoggleSwing.sliderBox("HUD Brightness", hudBrightnessSliderListener,
                (int) (settings.getPostprocessingHudBrightnessMin()),
                (int) (settings.getPostprocessingHudBrightnessMax()), 1,
                (int) (settings.getPostprocessingHudBrightness()), new JLabel("")));

        final ChangeListener sphereBrightnessSliderListener = new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                final JSlider source = (JSlider) e.getSource();
                if (source.hasFocus()) {
                    settings.setPostprocessingSphereBrightness(source.getValue());
                }
            }
        };
        targetPanel.add(GoggleSwing.sliderBox("Sphere Brightness", sphereBrightnessSliderListener,
                (int) (settings.getPostprocessingSphereBrightnessMin()),
                (int) (settings.getPostprocessingSphereBrightnessMax()), 1,
                (int) (settings.getPostprocessingSphereBrightness()), new JLabel("")));

        targetPanel.add(GoggleSwing.verticalStrut(1));
    }

    private void createPointCloudPanel(JPanel targetPanel) {
        final ItemListener pointGasDistanceDependantSizeListener = new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent e) {
                final JCheckBox source = (JCheckBox) e.getSource();
                settings.setPointGasPointSizeDependantOnCameraDistance(source.isSelected());
            }
        };

        final ChangeListener pointGasSizeSliderListener = new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                final JSlider source = (JSlider) e.getSource();
                if (source.hasFocus()) {
                    settings.setPointGasPointSizeSetting(source.getValue());
                }
            }
        };

        final ChangeListener pointGasBrightnessSliderListener = new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                final JSlider source = (JSlider) e.getSource();
                if (source.hasFocus()) {
                    settings.setPostprocessingPointGasBrightness(source.getValue());
                }
            }
        };

        final ChangeListener blurPassListener = new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                final JSlider source = (JSlider) e.getSource();
                if (source.hasFocus()) {
                    settings.setPointGasBlurPassSetting(source.getValue());
                }
            }
        };

        final ChangeListener blurTypeListener = new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                final JSlider source = (JSlider) e.getSource();
                if (source.hasFocus()) {
                    settings.setPointGasBlurTypeSetting(source.getValue());
                }
            }
        };

        final ChangeListener blurSizeListener = new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                final JSlider source = (JSlider) e.getSource();
                if (source.hasFocus()) {
                    settings.setPointGasBlurSizeSetting(source.getValue());
                }
            }
        };

        ActionListener scientificPresetListener = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                settings.setPointGasPointSizeDependantOnCameraDistance(false);
                settings.setPointGasPointSizeSetting(1);
                settings.setPointGasBlurPassSetting(1);
                settings.setPointGasBlurSizeSetting(2);
                settings.setPointGasBlurTypeSetting(1);
                settings.setPostprocessingPointGasBrightness(10);

                pointCloudConfig.removeAll();
                createPointCloudPanel(pointCloudConfig);
                validate();
                repaint();
            }
        };
        ActionListener embellishedPresetListener = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                settings.setPointGasPointSizeDependantOnCameraDistance(false);
                settings.setPointGasPointSizeSetting(1);
                settings.setPointGasBlurPassSetting(3);
                settings.setPointGasBlurSizeSetting(3);
                settings.setPointGasBlurTypeSetting(6);
                settings.setPostprocessingPointGasBrightness(25);

                pointCloudConfig.removeAll();
                createPointCloudPanel(pointCloudConfig);
                validate();
                repaint();

            }
        };
        ActionListener gaseousPresetListener = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                settings.setPointGasPointSizeDependantOnCameraDistance(false);
                settings.setPointGasPointSizeSetting(1);
                settings.setPointGasBlurPassSetting(6);
                settings.setPointGasBlurSizeSetting(7);
                settings.setPointGasBlurTypeSetting(8);
                settings.setPostprocessingPointGasBrightness(15);

                pointCloudConfig.removeAll();
                createPointCloudPanel(pointCloudConfig);
                validate();
                repaint();

            }
        };

        targetPanel.add(GoggleSwing.buttonBox("Presets", new GoggleSwing.ButtonBoxItem("Points",
                scientificPresetListener), new GoggleSwing.ButtonBoxItem("Blobs", embellishedPresetListener),
                new GoggleSwing.ButtonBoxItem("Gaseous", gaseousPresetListener)));

        targetPanel.add(GoggleSwing.checkboxBox("", new GoggleSwing.CheckBoxItem("Size Dependant on camera Distance",
                settings.isPointgasSizeDependantOnCameraDistance(), pointGasDistanceDependantSizeListener)));

        targetPanel.add(GoggleSwing.sliderBox("Point Gas Size", pointGasSizeSliderListener,
                (settings.getPointGasPointSizeMin()), (settings.getPointGasPointSizeMax()), 1,
                (settings.getPointGasPointSizeSetting()), new JLabel("")));
        targetPanel.add(GoggleSwing.sliderBox("Point Gas Brightness", pointGasBrightnessSliderListener,
                (int) (settings.getPostprocessingPointGasBrightnessMin()),
                (int) (settings.getPostprocessingPointGasBrightnessMax()), 1,
                (int) (settings.getPostprocessingPointGasBrightness()), new JLabel("")));
        targetPanel.add(GoggleSwing.sliderBox("Point Gas Blur Type", blurTypeListener, settings.getBlurTypeMin(),
                settings.getBlurTypeMax(), 1, settings.getPointGasBlurTypeSetting(), new JLabel("")));

        targetPanel.add(GoggleSwing.verticalStrut(1));
        targetPanel.add(GoggleSwing.sliderBox("Point Gas Blur Passes", blurPassListener, settings.getBlurPassMin(),
                settings.getBlurPassMax(), 1, settings.getPointGasBlurPassSetting(), new JLabel("")));

        targetPanel.add(GoggleSwing.verticalStrut(1));
        targetPanel.add(GoggleSwing.sliderBox("Point Gas Blur Size", blurSizeListener, settings.getBlurSizeMin(),
                settings.getBlurSizeMax(), 1, settings.getPointGasBlurSizeSetting(), new JLabel("")));
    }

    private void createStarPanel(JPanel targetPanel) {
        final ChangeListener particleSizeMultiplierListener = new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                final JSlider source = (JSlider) e.getSource();
                if (source.hasFocus()) {
                    settings.setParticleSizeMultiplier(source.getValue());
                }
            }
        };
        targetPanel.add(GoggleSwing.sliderBox("Star size multiplier.", particleSizeMultiplierListener,
                (settings.getParticleSizeMultiplierMin()), (settings.getParticleSizeMultiplierMax()), 1,
                (int) (settings.getParticleSizeMultiplier()), new JLabel("")));

        final ChangeListener starBrightnessSliderListener = new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                final JSlider source = (JSlider) e.getSource();
                if (source.hasFocus()) {
                    settings.setPostprocessingStarBrightness(source.getValue());
                }
            }
        };
        targetPanel.add(GoggleSwing.sliderBox("Star Brightness", starBrightnessSliderListener,
                (int) (settings.getPostprocessingStarBrightnessMin()),
                (int) (settings.getPostprocessingStarBrightnessMax()), 1,
                (int) (settings.getPostprocessingStarBrightness()), new JLabel("")));

        targetPanel.add(GoggleSwing.verticalStrut(1));

        final ChangeListener starHaloBrightnessSliderListener = new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                final JSlider source = (JSlider) e.getSource();
                if (source.hasFocus()) {
                    settings.setPostprocessingStarHaloBrightness(source.getValue());
                }
            }
        };
        targetPanel.add(GoggleSwing.sliderBox("Star Halo Brightness", starHaloBrightnessSliderListener,
                (int) (settings.getPostprocessingStarHaloBrightnessMin()),
                (int) (settings.getPostprocessingStarHaloBrightnessMax()), 1,
                (int) (settings.getPostprocessingStarHaloBrightness()), new JLabel("")));

        targetPanel.add(GoggleSwing.verticalStrut(1));

        final ChangeListener blurTypeListener = new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                final JSlider source = (JSlider) e.getSource();
                if (source.hasFocus()) {
                    settings.setStarHaloBlurTypeSetting(source.getValue());
                }
            }
        };
        targetPanel.add(GoggleSwing.sliderBox("Star Halo Blur Type", blurTypeListener, settings.getBlurTypeMin(),
                settings.getBlurTypeMax(), 1, settings.getStarHaloBlurTypeSetting(), new JLabel("")));

        targetPanel.add(GoggleSwing.verticalStrut(1));

        final ChangeListener blurPassListener = new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                final JSlider source = (JSlider) e.getSource();
                if (source.hasFocus()) {
                    settings.setStarHaloBlurPassSetting(source.getValue());
                }
            }
        };
        targetPanel.add(GoggleSwing.sliderBox("Star Halo Blur Passes", blurPassListener, settings.getBlurPassMin(),
                settings.getBlurPassMax(), 1, settings.getStarHaloBlurPassSetting(), new JLabel("")));

        targetPanel.add(GoggleSwing.verticalStrut(1));

        final ChangeListener blurSizeListener = new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                final JSlider source = (JSlider) e.getSource();
                if (source.hasFocus()) {
                    settings.setStarHaloBlurSizeSetting(source.getValue());
                }
            }
        };
        targetPanel.add(GoggleSwing.sliderBox("Star Halo Blur Size", blurSizeListener, settings.getBlurSizeMin(),
                settings.getBlurSizeMax(), 1, settings.getStarHaloBlurSizeSetting(), new JLabel("")));

    }

    private JPanel createBottomPanel() {
        final JPanel bottomPanel = new JPanel();
        bottomPanel.setFocusCycleRoot(true);
        bottomPanel.setFocusTraversalPolicy(new FocusTraversalPolicy() {
            @Override
            public Component getComponentAfter(Container aContainer, Component aComponent) {
                return null;
            }

            @Override
            public Component getComponentBefore(Container aContainer, Component aComponent) {
                return null;
            }

            @Override
            public Component getDefaultComponent(Container aContainer) {
                return null;
            }

            @Override
            public Component getFirstComponent(Container aContainer) {
                return null;
            }

            // No focus traversal here, as it makes stuff go bad (some things
            // react on focus).
            @Override
            public Component getLastComponent(Container aContainer) {
                return null;
            }
        });

        final JButton oneForwardButton = GoggleSwing.createImageButton("images/media-playback-oneforward.png", "Next",
                null);
        final JButton oneBackButton = GoggleSwing.createImageButton("images/media-playback-onebackward.png",
                "Previous", null);
        final JButton rewindButton = GoggleSwing.createImageButton("images/media-skip-backward.png", "Rewind", null);
        final JButton screenshotButton = GoggleSwing.createImageButton("images/camera.png", "Screenshot", null);
        final JButton playButton = GoggleSwing.createImageButton("images/media-playback-start.png", "Start", null);
        final ImageIcon playIcon = GoggleSwing.createImageIcon("images/media-playback-start.png", "Start");
        final ImageIcon stopIcon = GoggleSwing.createImageIcon("images/media-playback-stop.png", "Start");

        bottomPanel.setLayout(new BoxLayout(bottomPanel, BoxLayout.Y_AXIS));
        final JPanel buttonPanel = new JPanel();
        final JPanel timebarPanel = new JPanel();
        buttonPanel.setLayout(new BoxLayout(buttonPanel, BoxLayout.X_AXIS));
        timebarPanel.setLayout(new BoxLayout(timebarPanel, BoxLayout.X_AXIS));

        screenshotButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                timer.setScreenshotNeeded(true);
            }
        });
        buttonPanel.add(screenshotButton);

        buttonPanel.add(GoggleSwing.horizontalStrut(2));

        rewindButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                timer.rewind();
                playButton.setIcon(playIcon);
            }
        });
        buttonPanel.add(rewindButton);

        buttonPanel.add(GoggleSwing.horizontalStrut(2));

        oneBackButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                timer.oneBack();
                playButton.setIcon(playIcon);
            }
        });
        buttonPanel.add(oneBackButton);

        buttonPanel.add(GoggleSwing.horizontalStrut(2));

        playButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (timer.isPlaying()) {
                    timer.stop();
                    playButton.setIcon(playIcon);
                } else {
                    timer.start();
                    playButton.setIcon(stopIcon);
                }
            }
        });
        buttonPanel.add(playButton);

        buttonPanel.add(GoggleSwing.horizontalStrut(2));

        oneForwardButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                timer.oneForward();
                playButton.setIcon(playIcon);
            }
        });
        buttonPanel.add(oneForwardButton);

        buttonPanel.add(GoggleSwing.horizontalStrut(20));

        frameCounter = new JFormattedTextField();
        frameCounter.setValue(new Integer(1));
        frameCounter.setColumns(4);
        frameCounter.setMaximumSize(new Dimension(40, 20));
        frameCounter.setValue(0);
        frameCounter.addPropertyChangeListener(new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent e) {
                final JFormattedTextField source = (JFormattedTextField) e.getSource();
                if (source.hasFocus()) {
                    if (source == frameCounter) {
                        if (timer.isInitialized()) {
                            timer.setFrame(((Number) frameCounter.getValue()).intValue() - timeBar.getMinimum(), false);
                        }
                        playButton.setIcon(playIcon);
                    }
                }
            }
        });

        buttonPanel.add(frameCounter);

        timeBar.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                final JSlider source = (JSlider) e.getSource();
                if (source.isFocusOwner() && !source.getValueIsAdjusting()) {
                    timer.setFrame(timeBar.getValue(), false);
                    playButton.setIcon(playIcon);
                }
            }
        });
        timebarPanel.add(timeBar);

        bottomPanel.add(buttonPanel);
        bottomPanel.add(timebarPanel);

        return bottomPanel;
    }

    private void makeTimer() {
        if (timer.isInitialized()) {
            timer.close();
        }
        timer = new GlueTimedPlayer(timeBar, frameCounter);

    }

    public static GlueTimedPlayer getTimer() {
        return timer;
    }
}
