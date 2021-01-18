# Cancer Bioinformatics

Course material for the BIOS691 "Cancer Bioinformatics" course, January 25 - May 7, 2021. Course web site is https://bios691-cancer-bioinformatics.netlify.app/. Links to images and material are added to the slides. 

## Course Description

Title: "Cancer Bioinformatics"   
Course web site: https://bios691-cancer-bioinformatics.netlify.app/   
Time: Monday, Wednesday September 21 - October 23, 2020, 9:00am-10:20am    
Location: online only, Zoom (link will be provided)   
Equipment Required:  A computer with Internet access  
Registrant Limit: 25    
Description: This is a graduate course in cancer bioinformatics. This course is designed to provide an introduction to genomics, sequencing technologies, data analysis. We will cover the technology and statistics of omics data, common statistical methods and tools for data analysis and visualization. Knowledge of basic statistics and R programming is a plus but not required.  

[Other Computational Genomics Courses](https://bios691-cancer-bioinformatics.netlify.app/reading/)

## Previous courses

[Statistical Methods for High-throughput Genomic Data I](https://mdozmorov.github.io/BIOS567.2018/)  
[Statistical Methods for High-throughput Genomic Data II](https://mdozmorov.github.io/BIOS668.2018/)

## Site template

Site template is based on the ["Statistical Image Analysis Course for Neuroscientists"](https://github.com/laderast/stats_for_neuroscientists) repository developed by [Ted Laderas](https://laderast.github.io/). The original source is ["GSU MPA/MPP course on program evaluation and causal inference"](https://github.com/andrewheiss/evalsp20.classes.andrewheiss.com) repository by [Andrew Heiss](https://www.andrewheiss.com/)

## Site building 

```
# Build and serve site
blogdown::serve_site()
# Stop server
blogdown::stop_server()
```

## Theme

This site uses the [Academic Hugo theme](https://sourcethemes.com/academic/), with some slight template modifications found in `/assets/` and `layouts/`. The theme is included as a submodule, so when when cloning for the first time, use this command to get the theme too:

    git clone --recursive https://github.com/gcushen/hugo-academic.git

To get the theme later, use this command:

    git submodule add \
      https://github.com/gcushen/hugo-academic.git \
      themes/hugo-academic

To update to the latest version of the theme, use:

    git submodule update --recursive --remote
