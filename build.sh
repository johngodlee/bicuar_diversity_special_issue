#!/bin/bash
{

IMG="manuscript/img/"
INC="manuscript/include/"

# Build script for manuscript by running each file in the right order

# Remove intermediate objects
# rm data/trees_ab_mat.rds
# rm data/trees_list.rds
# rm data/stems_list.rds
# rm data/plot_clim_list.rds
# rm data/plot_div_list.rds
# rm data/plot_loc_list.rds
# rm data/region_temp_precip.rds

# Run data analysis
#Rscript scripts/temp_precip_region.R
#Rscript scripts/stems_clean.R
#Rscript scripts/species_clean.R
#Rscript scripts/bicuar_20_50.R
#Rscript scripts/stems_to_trees.R
#Rscript scripts/plot_loc.R
#Rscript scripts/plot_clim.R
#Rscript scripts/plot_div.R
#Rscript scripts/beta_div.R
#Rscript scripts/data_descrip.R
#Rscript scripts/dbh_bin.R
#Rscript scripts/degrad.R
#Rscript scripts/bicuar_map.R
#Rscript scripts/bicuar_species_table.R
#Rscript scripts/plot_clim_plot.R

# Copy analysis outputs to manuscript directory
mkdir -p $IMG

cp img/plot_map.pdf $IMG
cp img/temp_precip.pdf $IMG
cp img/all_nmds_envfit.pdf $IMG
cp img/div_box.pdf $IMG
cp img/degrad_box.pdf $IMG
cp img/bicuar_degrad_nmds.pdf $IMG
cp img/bicuar_map.pdf $IMG
cp img/stem_ab_dbh_bin_group.pdf $IMG

mkdir -p $INC

cp include/group_descrip.tex $INC
cp include/anova_table.tex $INC
cp include/anova_degrad.tex $INC
cp include/site_pairs_js.tex $INC
cp include/bicuar_species.tex $INC
cp include/data_descrip_figures.tex $INC
cp include/dbh_bin_figures.tex $INC
cp include/beta_div_figures.tex $INC
cp include/plot_div_figures.tex $INC
cp include/degrad_figures.tex $INC

# Edit tex inputs 
# \& {$H'\''$} \& {Basal area (m\\textsuperscript{2} ha\\textsuperscript{-1})}
sed -i 's/\\extracolsep{5pt}} cccccc/\\extracolsep{0pt}} rccccc/g' "${INC}group_descrip.tex"
sed -i '10s/.*/{Plot group} \& \\multicolumn{1}{p{1cm}}{\\centering MAT \\\\ ($^\\circ$C)} \& \\multicolumn{1}{p{1.5cm}}{\\centering MAP \\\\ (mm y\\textsuperscript{-1})} \& \\multicolumn{1}{p{1.5cm}}{\\centering CWD \\\\ (mm y\\textsuperscript{-1})} \& \\multicolumn{1}{p{2cm}}{\\centering Latitude \\\\ (DD)} \& \\multicolumn{1}{p{2cm}}{\\centering Longitude \\\\ (DD)} \& {N plots} \& {N species} \\\\/' "${INC}group_descrip.tex"
sed -i 's/caption{}/caption{Description of each group of plots used in the analysis. MAT = Mean Annual Temperature, MAP = Mean Annual Precipitation, CWD = Climatic Water Deficit, DD = Decimal Degrees.}/g' "${INC}group_descrip.tex"

sed -i 's/bicuar/Angola (ANG)/g' "${INC}group_descrip.tex"
sed -i 's/drc/DRC (DRC)/g' "${INC}group_descrip.tex"
sed -i 's/kilwa/Tanzania (TZA)/g' "${INC}group_descrip.tex"
sed -i 's/nham/Mozambique (MOZ)/g' "${INC}group_descrip.tex"

sed -i 's/\\extracolsep{5pt} ccccc/\\extracolsep{0pt} lcccc/g' "${INC}nmds_envfit.tex"
sed -i '10s/.*/{Var.} \& {NMDS1} \& {NMDS2} \& {R\\textsuperscript{2}} \& {Prob.} \\\\/' "${INC}nmds_envfit.tex"
sed -i 's/caption{}/caption{Test statistics for environmental fits on the NMDS of plot species composition.}/g' "${INC}nmds_envfit.tex"
sed -i '15s/0\.7/0.70/' "${INC}nmds_envfit.tex"

sed -i 's/\\extracolsep{5pt}} cccc/\\extracolsep{0pt}} rrcc/g' "${INC}site_pairs_js.tex"
sed -i 's/bicuar/Angola/g' "${INC}site_pairs_js.tex"
sed -i 's/drc/DRC/g' "${INC}site_pairs_js.tex"
sed -i 's/kilwa/Tanzania/g' "${INC}site_pairs_js.tex"
sed -i 's/nham/Mozambique/g' "${INC}site_pairs_js.tex"
sed -i '10s/.*/Site 1 \& Site 2 \& $S_{S}$ \& Shared species \\\\/g' "${INC}site_pairs_js.tex"
sed -i 's/caption{}/caption{Pairwise beta diversity comparison of plot groups measured by the S\\o{}rensen coefficient ($S_s$) of percentage similarity of aggregated plot level data from each of the four sites. Values in brackets are the number of species unique to each site in each comparison.}/g' "${INC}site_pairs_js.tex"

sed -i 's/textbackslash //g' "${INC}bicuar_species.tex"
sed -i 's/\\{/{/g' "${INC}bicuar_species.tex"
sed -i 's/\\}/}/g' "${INC}bicuar_species.tex"
sed -i '10s/.*/{Family} \& {Species} \& \\multicolumn{1}{p{2cm}}{\\centering Stem diam. \\\\ (cm)} \& \\multicolumn{1}{p{2cm}}{\\centering Basal area \\\\ (m\\textsuperscript{2} ha\\textsuperscript{-1})} \& {N stems} \& {N stems ha\\textsuperscript{-1}} \\\\/' "${INC}bicuar_species.tex"
sed -i 's/\\extracolsep{5pt}} cccccc/\\extracolsep{-5pt}} rrcccc/g' "${INC}bicuar_species.tex"
sed -i 's/caption{}/caption{Species found in one hectare plots in Bicuar National Park. Stem diameter and basal area are the mean of all stems with the standard error of the mean in brackets. Number of stems per hectare is mean of the number of stems in all one hectare plots where stems of that species are present with the standard error of the mean in brakcets. Species found only in Bicuar National Park are marked in bold text with an asterisk.}/g' "${INC}bicuar_species.tex"

sed -i 's/\\extracolsep{5pt}}lcccc/\\extracolsep{0pt}}lcccc/g' "${INC}anova_table.tex"
sed -i 's/groupdrc/DRC/'"${INC}anova_table.tex"
sed -i 's/groupkilwa/Tanzania/' "${INC}anova_table.tex"
sed -i 's/groupnham/Mozambique/' "${INC}anova_table.tex"
sed -i 's/rich/Species richness/' "${INC}anova_table.tex"
sed -i 's/basal\\_area/Basal area/g' "${INC}anova_table.tex"
sed -i "s/shannon/Shannon ($H'$)/g" "${INC}anova_table.tex"
sed -i 's/shannon\\_equit/Shannon equitability ($E_{H}$)/g' "${INC}anova_table.tex"
sed -i '6d' "${INC}anova_table.tex"
sed -i 's/caption{}/caption{Results of ANOVA tests for alpha diversity metrics and plot basal area, among the four sites. Mean values for each site with standard errors in parentheses are shown. Asterisks indicate the p-value of individual sites (***<0.001, **<0.01, *<0.05, .<0.1).}/g' "${INC}anova_table.tex"
sed -i '33d' "${INC}anova_table.tex"


# Compile manuscript
cd manuscript
# pdflatex special_issue_manuscript.tex
# pdflatex special_issue_manuscript.tex
# bibtex special_issue_manuscript
# bibtex special_issue_manuscript
# pdflatex special_issue_manuscript.tex
latexmk -C special_issue_manuscript.tex

# Compile README.md
cd ..
paper_title="$(grep "Title" manuscript/special_issue_manuscript.tex | awk -F "{|}" '{print $2}')"
abstract="$(grep -e "abstract{" manuscript/special_issue_manuscript.tex | awk -F "{|}" "{print $2}" | sed 's/\\abstract{//g' | sed 's/}$//g' )"
url="PLACEHOLDER"

cat >README.md <<EOL 
Code for "[${paper_title}](${url})"

Godlee, J. L.; Gonçalves, F. M.; Tchamba, J. J.; Chisingui, A. V.; Muledi, J. I.; Shutcha, M. N.; Ryan, C. M.; Dexter, K. G. ${paper_title}. _Diversity_ __2020__, 12(1), 1-19. 

> "${abstract}"

\`\`\`
@article{Godlee2020diversity,
  title = {${paper_title}},
  author = {Godlee, J. L. and Gon\\c{c}alves, F. M. and Tchamba, J. J. and Chisingui, A. V. and Muledi, J. I. and Shutcha, M. N. and Ryan, C. M. and Dexter, K. G.},
  journal = {Diversity},
  year = 2020,
  volume = {12},
  issue = {1},
  pages = {}
}
\`\`\`

* \`manuscript/\` contains the LaTeX build files for the paper
* \`data/\` contains the data used in analyses
* \`scripts/\` contains R scripts used for analysis
* \`img/\` contains all images generated by the R scripts
* \`include/\` contains all non-image outputs generated by the R scripts
* \`cover_letter/\` contains the cover letter sent to guest editors on submission
* \`misc/\` contains a number of obsolete files from previous versions of the manuscript 
* \`build.sh\` is a shell script which runs analysis and builds the paper in the corect order

EOL




} > output.log 2>&1

