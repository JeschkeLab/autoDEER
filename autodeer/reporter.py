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
            template = PageTemplate(id='normal', frames=frame, onPage=self.header, onPageEnd=self.footer)
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
