# -*- coding: utf-8 -*-
"""
    sphinx.directives.code
    ~~~~~~~~~~~~~~~~~~~~~~

    :copyright: Copyright 2007-2011 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import sys
import codecs

from docutils.nodes import Body, Element
from docutils import nodes
#from docutils.parsers.rst import Directive, directives
from docutils.parsers.rst import directives

from sphinx import addnodes
from sphinx.util import parselinenos
from sphinx.util.nodes import set_source_info
from sphinx.util.compat import Directive


class sagecell(Body, Element):
   pass

class SageCell(Directive):
    """
    Directive for a code block with special highlighting or line numbering
    settings.
    """

    has_content = True
    required_arguments = 0
    optional_arguments = 0
    final_argument_whitespace = False
    option_spec = {
        'linenos': directives.flag,
        'emphasize-lines': directives.unchanged_required,
    }

    def run(self):
        code = u'\n'.join(self.content)

        linespec = self.options.get('emphasize-lines')
        if linespec:
            try:
                nlines = len(self.content)
                hl_lines = [x+1 for x in parselinenos(linespec, nlines)]
            except ValueError, err:
                document = self.state.document
                return [document.reporter.warning(str(err), line=self.lineno)]
        else:
            hl_lines = None

        literal = nodes.literal_block(code, code)
        literal['language'] = 'python'
        literal['linenos'] = 'linenos' in self.options
        if hl_lines is not None:
            literal['highlight_args'] = {'hl_lines': hl_lines}
        set_source_info(self, literal)
        return [literal]


def setup(app):
    app.add_directive('sagecell', SageCell)


#directives.register_directive('sagecell', SageCell)
