from reportlab.platypus import SimpleDocTemplate, Flowable, Preformatted, Spacer, PageBreak
from reportlab.graphics.shapes import Drawing
from reportlab.lib import pagesizes
from reportlab.lib.units import cm
from reportlab.platypus.paragraph import Paragraph
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.graphics import renderPDF
from reportlab.platypus.doctemplate import PageTemplate
from reportlab.platypus.frames import Frame
from functools import partial
from datetime import date
import matplotlib.pyplot as plt
from autodeer.DEER_analysis import DEERanalysis_plot, plot_overlap

from svglib.svglib import svg2rlg
from io import BytesIO
import os

from collections import OrderedDict
styles = getSampleStyleSheet()


package_dir = os.path.dirname(os.path.split(__file__)[0])

class SvgFlowable(Flowable):
    """Convert byte streasm containing SVG into a Reportlab Flowable."""

    def __init__(self, svg: BytesIO) -> None:
        """Convert SVG to RML drawing on initializtion."""
        svg.seek(0)
        self.drawing: Drawing = svg2rlg(svg)
        self.width: int = self.drawing.minWidth()
        self.height: int = self.drawing.height
        self.drawing.setProperties({"vAlign": "CENTER", "hAlign": "CENTER"})

    def wrap(self, *_args):
        """Return diagram size."""
        return (self.width, self.height)

    def draw(self) -> None:
        """Render the chart."""
        renderPDF.draw(self.drawing, self.canv, 0, 0)


class Reporter():

    def __init__(self, filepath, pagesize='A4') -> None:

            # Build the pdf

            if pagesize == 'A4':
                PAGESIZE = pagesizes.portrait(pagesizes.A4)
            elif pagesize =='Letter':
                PAGESIZE = pagesizes.portrait(pagesizes.LETTER)
            else:
                raise ValueError("Only pagesizes of 'A4' or 'Letter' are supported")
            
            self.pdf = SimpleDocTemplate(filepath, pagesize=PAGESIZE, 
                    leftMargin = 2.2 * cm, 
                    rightMargin = 2.2 * cm,
                    topMargin = 2.5 * cm, 
                    bottomMargin = 2.5 * cm)
            
            self.story = OrderedDict()  # possibly change to a normal dict in the future
            frame = Frame(self.pdf.leftMargin, self.pdf.bottomMargin, self.pdf.width, self.pdf.height, id='normal')
            template = PageTemplate(id=None, frames=frame, onPage=self.header, onPageEnd=self.footer)
            self.pdf.addPageTemplates([template])

            pass
    
    def header(self, canvas, doc):
        logo = os.path.abspath(package_dir + '/docsrc/_static/autoDEER_light.png')
        canvas.saveState()
        canvas.setFont('Times-Bold',12)
        date_str = date.today().strftime("%d/%m/%Y")
        canvas.drawImage(logo, 10, doc.height+80, width=215, height=50,mask='auto')

        canvas.drawCentredString(doc.width-2, doc.height+95, f'autodeer report: {date_str}')
        canvas.restoreState()

    def footer(self, canvas, doc):
        canvas.saveState()
        page_number = canvas.getPageNumber()
        canvas.setFont('Times-Roman',9)
        canvas.drawString(2.2 * cm, 2.2 * cm, f'Page {page_number}')
        canvas.restoreState()

    


    def _build(self):
        # Convert ordered dict to list of values and then flatten

        flat_list = []
        for key in self.story:
            flat_list += self.story[key]

        self.pdf.build(flat_list)

    def add_title(self, key, title):
        self.story[key] = [Paragraph(title, styles['Title'])]

    def add_new_section(self, key, title):
        substory = []
        substory.append(Paragraph(title, styles['Heading2']))
        self.story[key] = substory

    def add_text(self, key, text, title = None):
        if title is not None:
            self.story[key].append(Paragraph(title, styles['Heading4']))
        self.story[key].append(Paragraph(text, styles['Normal']))
    
    def add_code_block(self, key, code, title = None):
        if title is not None:
            self.story[key].append(Paragraph(title, styles['Heading4']))
        self.story[key].append(Preformatted(code, styles['Code']))
    
    def add_figure(self, key, figure, title = None):
        # Add a matplotlib figure to the reportlab document
        if title is not None:
            self.story[key].append(Paragraph(title, styles['Heading4']))

        imgdata = BytesIO()
        # Resize figure whilst preservaing aspect ratio
        width = figure.get_figwidth()
        height = figure.get_figheight()
        aspect = height/width
        New_width = 5
        figure.set_figwidth(New_width)
        figure.set_figheight(New_width*aspect)
        figure.tight_layout()
        figure.savefig(imgdata, format='svg')
        imgdata.seek(0)  # rewind the data

        self.story[key].append(SvgFlowable(imgdata))

    def add_space(self, key, height=5):
        self.story[key].append(Spacer(width=100, height=height))

    
    def add_page_break(self, key):
        self.story[key].append(PageBreak())

    def add_table(self, key, lists):
        """Generates a table as a reportlab flowable from a list of lists
        """
        t = Table(lists)
        t.setStyle(TableStyle([('ALIGN',(1,1),(-2,-2),'RIGHT'),]))
        self.story[key].append(t)   


def create_report(save_path, Results:dict, SpectrometerInfo:dict=None, UserInputs:dict=None, Pulses=None):
        report = Reporter(filepath=save_path,pagesize='A4')

        report.add_title('title','autoDEER Report')
        if SpectrometerInfo is not None:
            report.add_new_section('spec',' Spectrometer')
            report.add_text('spec', 'Local Name: ' + SpectrometerInfo['Local Name'])
            report.add_text('spec', 'Manufacturer: ' + SpectrometerInfo['Manufacturer'])
            report.add_text('spec', 'Model: ' + SpectrometerInfo['Model'])

            report.add_text('spec', 'Resonator: ' + SpectrometerInfo['Resonator'])
            report.add_space('spec', height=10)
        report.add_new_section('inputs',' User Inputs')
        if UserInputs is not None:
            report.add_text('inputs', 'Sample Name: ' + UserInputs['Sample Name'])
            report.add_text('inputs', f"Temperature: {UserInputs['Temp']:.3f} K")
            report.add_text('inputs', f"Max Time: {UserInputs['MaxTime']} hrs")

        report.add_page_break('inputs')

        if 'fieldsweep' in Results:
            report.add_new_section('fieldsweep',' Field Sweep')
            fig,axs = plt.subplots(1,1,figsize=(5, 3))
            Results['fieldsweep'].plot(axs=axs, fig=fig)
            report.add_figure('fieldsweep', fig)
            report.add_text('fieldsweep', f"Gyromagnetic Ratio: {Results['fieldsweep'].gyro:.3g} GHz/G")
            if hasattr(Results['fieldsweep'], 'results'):
                report.add_space('fieldsweep')
                report.add_code_block('fieldsweep', Results['fieldsweep'].results.__str__(), title='Fit Results')
            report.add_page_break('fieldsweep')

        if 'respro' in Results:
            report.add_new_section('respro',' Resonator Profile')
            fig,axs = plt.subplots(1,1,figsize=(5, 5))
            fitresult = Results['respro']
            if 'fieldsweep'in Results:
                fitresult.plot(fieldsweep=Results['fieldsweep'],axs=axs, fig=fig);
            else:
                fitresult.plot(axs=axs,fig=fig)

            report.add_figure('respro', fig)
            if hasattr(Results['respro'], 'results'):
                report.add_space('respro')
                report.add_code_block('respro', Results['respro'].results.__str__(), title='Fit Results')
            report.add_page_break('respro')

        if Pulses is not None:
            pump_pulse = Pulses['pump_pulse']
            exc_pulse = Pulses['exc_pulse']
            ref_pulse = Pulses['ref_pulse']
            report.add_new_section('pulses',' Optimised Pulses')
            fig,axs = plt.subplots(1,1,figsize=(5, 5))
            plot_overlap(Results['fieldsweep'], pump_pulse, exc_pulse,ref_pulse, axs=axs,fig=fig)

            report.add_figure('pulses', fig)
            report.add_space('pulses')

            report.add_code_block('pulses', exc_pulse.__str__(), title='Excitation Pulse')
            report.add_code_block('pulses', ref_pulse.__str__(), title='Refocusing Pulse')
            report.add_code_block('pulses', pump_pulse.__str__(), title='Pump Pulse')


        if 'relax' in Results:
            report.add_new_section('relax',' Relaxation')
            fig,axs = plt.subplots(1,1,figsize=(5, 5))
            Results['relax'].plot(axs=axs, fig=fig)
            report.add_figure('relax', fig)
            if hasattr(Results['relax'], 'results'):
                report.add_space('relax')
                report.add_code_block('relax', Results['relax'].results.__str__(), title='Fit Results')
            report.add_page_break('relax') 


        if 'quickdeer' in Results:
            report.add_new_section('quickdeer',' QuickDEER')
            fig,axs = plt.subplot_mosaic([['Primary_time'], 
                                          ['Primary_dist']], figsize=(6,6))
            
            DEERanalysis_plot(Results['quickdeer'], background=True, ROI=Results['quickdeer'].ROI, axs= axs, fig=fig,text=False)
            report.add_figure('quickdeer', fig)
            
            if 'quickdeer' in Results:
                report.add_space('quickdeer')
                report.add_code_block('quickdeer', Results['quickdeer'].__str__(), title='Fit Results')
            report.add_page_break('quickdeer')
            
        report._build()
        pass

def combo_figure(EDFS, respro, pulses:dict, relaxation:list, init_deer, long_deer ,title=None):
    """
    Creates a 2x2 summary figure. 
        - The top left plot is the EDFS and resonator profile, overlapped with the optimised pulses. 
        - The top right plot is the relaxation data and fits.
        - The bottom left plot is the initial DEER data and fits.
        - The bottom right plot is the final DEER data and fits.
    
    Parameters
    ----------
    EDFS: ad.FieldSweepAnalysis
        The Echo-Detected Field Sweep analysis.
    respro: ad.ResonatorProfileAnalysis
        The resonator profile analysis
    pulses: dict
        A dictionary containing the optimised pulses.
    relaxation: list
        A list containing the relaxation data and fits.
    init_deer: deerlab.FitResult
        The initial DEER data and fits.
    long_deer: deerlab.FitResult
        The final DEER data and fits.
    title: str, optional
        The title of the figure, by default None
    
    
    """
    figure = plt.figure(figsize=(10, 10),constrained_layout=True)
    figs = figure.subfigures(2, 2, height_ratios=(4,6), width_ratios=(1,1),wspace=.12)

    if title is not None:
        figure.suptitle(title,size=20)
    
    figs[0,0].suptitle('a. EDFS',size=15)
    plot_overlap(EDFS, pulses['pump_pulse'], pulses['exc_pulse'],pulses['ref_pulse'],respro=respro,fig=figs[0,0]);
    figs[0,0].axes[0].set_xlim(-0.3,0.1)
    figs[0,1].suptitle('b. Relaxation',size=15)
    plot_1Drelax(*relaxation,figs[0,1])
    figs[1,0].suptitle('c. Intial DEER',size=15)
    DEERanalysis_plot_pub(init_deer[0],figs[1,0]);
    figs[1,1].suptitle('d. Final DEER',size=15)
    DEERanalysis_plot_pub(long_deer,figs[1,1]);

    return figure