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
import nl.esciencecenter.asterisk.interfaces.TimedPlayer;
import nl.esciencecenter.esight.swing.CustomJSlider;
import nl.esciencecenter.esight.swing.GoggleSwing;

public class AsteriskInterfaceWindow extends JPanel {
    public static enum TweakState {
        NONE, DATA, VISUAL
    }

    private final AsteriskSettings settings = AsteriskSettings.getInstance();

    private static final long serialVersionUID = 1L;

    protected CustomJSlider timeBar;

    protected JFormattedTextField frameCounter;

    public static TimedPlayer timer;

    private final JTabbedPane configPanel;

    private final JPanel starConfig, pointCloudConfig, miscConfig;

    public AsteriskInterfaceWindow() {
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

        configPanel.setVisible(true);

        add(bottomPanel, BorderLayout.SOUTH);

        // Read command line file information
        makeTimer();

        // setTweakState(TweakState.VISUAL);

        setVisible(true);
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

    public static TimedPlayer getTimer() {
        return timer;
    }
}
