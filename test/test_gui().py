import pytest
from autodeer.gui import autoDEERUI
from PyQt6.QtWidgets import QApplication, QMainWindow
from autodeer.gui.autoDEER_worker import autoDEERWorker
from autodeer.hardware.dummy import dummyInterface
import time
from PyQt6 import QtCore, QtWidgets
import numpy as np



@pytest.fixture
def app(qtbot):
    test_app = QtWidgets.QApplication([])
    main = autoDEERUI()
    qtbot.addWidget(main)
    return main


def test_RunFullyAutoDEER(app, qtbot):
    # Set up user input values
    app.MaxTime.setValue(10)
    app.SampleName.setText('Test Sample')
    app.TempValue.setValue(25)

    # Connect to spectrometer
    app.load_spectrometer_config('config_files/Dummy.yaml')
    app.connect_spectrometer()

    app.spectromterInterface.speedup = 10000
    assert app.connected == True

    # Run FullyAutoDEER method
    worker = app.RunFullyAutoDEER()
    Running = True
    def set_running_false():
        Running = False
    def test_fieldsweep():
        qtbot.waitUntil(lambda : 'fieldsweep' in app.current_results, timeout= 240*1e3)
        time.sleep(1)

        assert app.gyroSpinBox.value() == pytest.approx(2.81079,rel=1e-2)
        assert app.gxSpinBox.value() == pytest.approx(2.0083,rel=1e-2)
        assert app.gySpinBox.value() == pytest.approx(2.0061,rel=1e-2)
        assert app.gzSpinBox.value() == pytest.approx(2.0021,rel=1e-2)
        assert app.AxSpinBox.value() == pytest.approx(13.68,rel=1)
        assert app.AySpinBox.value() == pytest.approx(13.68,rel=1)
        assert app.AzSpinBox.value() == pytest.approx(3.66,rel=1e-1)
        assert app.GBSpinBox.value() == pytest.approx(0.45,rel=1e-1)
        assert app.BoffsetSpinBox.value() == pytest.approx(0,rel=1e-1)

        fig = app.fsweep_canvas.figure

        assert len(fig.axes[0].lines) == 3

    def test_respro():
        qtbot.waitUntil(lambda : 'respro' in app.current_results, timeout= 120*1e3)

        assert app.centreFrequencyDoubleSpinBox.value() == pytest.approx(34.0,rel=1e-1)
        assert app.qDoubleSpinBox.value() == pytest.approx(80,rel=2e1)

        fig = app.respro_canvas.figure

        assert len(fig.axes[0].lines) == 1

    def test_relax_and_pulses():
        # test pulses
        assert app.ExcTypeLine.text() == 'Rect'
        assert app.ExcFreqBox.value() == pytest.approx(0,rel=1)
        assert app.ExcBWBox.value() == pytest.approx(16,rel=1)
        assert app.ExcBWBox.suffix() == ' ns'

        assert app.RefTypeLine.text() == 'Rect'
        assert app.RefFreqBox.value() == pytest.approx(0,rel=1)
        assert app.RefBWBox.value() == pytest.approx(16,rel=1)
        assert app.RefBWBox.suffix() == ' ns'

        assert app.PumpTypeLine.text() == 'HS'
        assert app.PumpFreqBox.value() == pytest.approx(-136,rel=10)
        assert app.PumpBWBox.value() == pytest.approx(220,rel=10)
        assert app.PumpBWBox.suffix() == ' MHz'


        qtbot.waitUntil(lambda : 'relax' in app.current_results, timeout= 120*1e3)

        assert app.DipolarEvoMax.value() == pytest.approx(2.1,rel=3e-1)
        assert app.DipolarEvo2hrs.value() == pytest.approx(1.9,rel=3e-1)

        fig = app.relax_canvas.figure

        assert len(fig.axes[0].lines) == 2

    def test_quickdeer_result():
        qtbot.waitUntil(lambda : 'quickdeer' in app.current_results, timeout= 10*60*1e3)
        
        assert app.q_DEER.ExperimentcomboBox.currentText() == '5pDEER'
        assert app.q_DEER.Tau1doubleSpinBox.value() == pytest.approx(4.3,rel=0.5)
        assert app.q_DEER.Tau2doubleSpinBox.value() == pytest.approx(4.3,rel=0.5)
        assert app.q_DEER.Tau3doubleSpinBox.value() == pytest.approx(0.3,rel=0.1)
        assert app.q_DEER.CompactnessradioButton.isChecked() == True
        assert app.q_DEER.PulseLengthdoubleSpinBox.value() == pytest.approx(128,rel=1)

        fig = app.q_DEER.static_canvas.figure

        assert len(fig.axes[0].lines) == 2
        assert len(fig.axes[1].lines) == 1

        fitresult = app.q_DEER.fitresult

        # count number of lam* in fitresult
        count = 0
        for key in fitresult.params.keys():
            if key.startswith('lam'):
                count += 1
        assert count == 10 # Both normal and uncertainty

        assert 'rec_dt' in fitresult
        assert 'rec_tau_max' in fitresult

        assert fitresult.ROI[0] == pytest.approx(5,rel=0.5)
        assert fitresult.ROI[1] == pytest.approx(3.5,rel=0.5)
        
    worker.signals.finished.connect(lambda: set_running_false())
    worker.signals.fsweep_result.connect(lambda: test_fieldsweep())
    worker.signals.respro_result.connect(lambda: test_respro())
    worker.signals.relax_result.connect(lambda: test_relax_and_pulses())
    worker.signals.quickdeer_result.connect(lambda: test_quickdeer_result())

    qtbot.waitUntil(lambda : Running == False, timeout=25*60*1e3)

    # Check that FullyAutoButton and AdvancedAutoButton are enabled again
    assert app.FullyAutoButton.isEnabled() == True
    assert app.AdvancedAutoButton.isEnabled() == True

