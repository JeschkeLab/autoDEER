import os
from autodeer.reporter import Reporter, SvgFlowable
from reportlab.lib.pagesizes import A4, LETTER
from reportlab.lib.units import cm
from reportlab.platypus import SimpleDocTemplate, Paragraph, Preformatted, Spacer, PageBreak
from reportlab.lib.styles import getSampleStyleSheet
from datetime import date
import pytest

@pytest.fixture
def reporter():
    return Reporter('test_report.pdf')

def test_reporter_init(reporter):
    assert isinstance(reporter.pdf, SimpleDocTemplate)
    assert isinstance(reporter.story, dict)

def test_reporter_add_title(reporter):
    reporter.add_title('test', 'Test Title')
    assert isinstance(reporter.story['test'][0], Paragraph)

def test_reporter_add_new_section(reporter):
    reporter.add_new_section('test', 'Test Section')
    assert isinstance(reporter.story['test'][0], Paragraph)

def test_reporter_add_text(reporter):
    reporter.add_new_section('test', 'Test Section')
    reporter.add_text('test', 'Test Text')
    assert isinstance(reporter.story['test'][0], Paragraph)

def test_reporter_add_code_block(reporter):
    reporter.add_new_section('test', 'Test Section')
    reporter.add_code_block('test', 'print("Hello, world!")')
    assert isinstance(reporter.story['test'][1], Preformatted)

def test_reporter_add_figure(reporter):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot([1, 2, 3], [4, 5, 6])
    reporter.add_new_section('test', 'Test Section')
    reporter.add_figure('test', fig)
    assert isinstance(reporter.story['test'][1], SvgFlowable)

def test_reporter_add_space(reporter):
    reporter.add_new_section('test', 'Test Section')
    reporter.add_space('test', height=10)
    assert isinstance(reporter.story['test'][1], Spacer)

def test_reporter_add_page_break(reporter):
    reporter.add_new_section('test', 'Test Section')
    reporter.add_page_break('test')
    assert isinstance(reporter.story['test'][1], PageBreak)

def test_reporter_build(reporter):
    reporter.add_new_section('test', 'Test Section')
    reporter.add_title('test', 'Test Title')
    reporter.add_new_section('test', 'Test Section')
    reporter.add_text('test', 'Test Text')
    reporter.add_code_block('test', 'print("Hello, world!")')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot([1, 2, 3], [4, 5, 6])
    reporter.add_figure('test', fig)
    reporter.add_space('test', height=10)
    reporter.add_page_break('test')
    reporter._build()
    assert os.path.exists('test_report.pdf')
    os.remove('test_report.pdf')