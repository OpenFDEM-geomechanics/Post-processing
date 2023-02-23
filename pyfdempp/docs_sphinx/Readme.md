# About Documentation
Documentation was put together using `Sphinx`, pointing at the .py files in openfdem folder. 

# Where to Find Config Files
All build settings for documentation are controlled in [conf.py](conf.py). The start page and display menu are found under the file [index.rst](index.rst).

# Running Sphinx

## HTML Build
To rerun the HTML build (ex. if a docstring is updated), run
 `make html` 
 in terminal.

To view the `HTML` files locally, navigate to  `_build` -> `html` -> [index.html](_build/html/index.html). From there the sub-modules can be accessed in your browser.
To format the HTML files, refer to `index.rst` and `conf.py`.

## Latex and PDF Build
To rerun the Latex and PDF build (ex. if a docstring is updated), run
 `make latexpdf` 
 in terminal.

If you encounter [Makefile:20: pdflatex] Error 2, run
`sphinx-build -M latexpdf . _build`


To view the PDF file, navigate to  `_build` -> `latex` -> [openfdempost-processing.pdf](_build/latex/openfdempost-processing.pdf). 

## pybadges

To insert badges on the top of the README files using the [pybadges pacakge](https://github.com/google/pybadges).

```console
python3 -m pybadges\
     --left-text="python"\
     --right-text="3.5|3.6|3.7|3.8|3.9|3.10"\
     --browser \
     --logo='https://dev.w3.org/SVG/tools/svgweb/samples/svg-files/python.svg'\
     --embed-logo yes
```

