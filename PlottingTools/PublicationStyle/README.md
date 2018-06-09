# Restyle Varial Figures for Publication

1. Install the `cms_figure` package.

   ```bash
   ./setup
   ```

2. Run the `restyle_2016.py` or `restyle_2017.py` scripts. The positional command line
   arguments are the path to the Varial output directory, the path to the new output
   directory for the restyled plots, and an optional `--logy` flag to produce restyled
   plots with a logarithmically scaled y-axis.

   For example,
   ```python
   python restyle_2016.py ../PlotWithVarial/my_varial_plots my_restyled_plots
   ```

3. I usually organize my plots into subdirectories within PlotWithVarial and the run the
   scripts over them. As an example, I've collected the commands I use into the bash script
   `style_plots.sh`.

