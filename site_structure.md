
- `config`
    - `_default`
        - `menus.yaml` - order and weights of menu items in the header (main) and the left side (assignment) panels
        - `params.yaml` - personalization settings

- `content`
    - `assignment` - content to appear under "Assignments" main menu
        - `_index.Rmd` - main assignment instructions. Preparation, problem sets, evaluation, code-through, exams, final project
        - `01-problem-set.Rmd` - description of the first problem set

    - `authors`
        - `MGD` - name initials
            - `_index.md` - about the author paragraph
    
    - `class` - content to appear under "Class" main menu
        - `_index.Rmd` - general class details description
        - `01-class.Rmd` - individual class/lab description, with embedded slides
    
    - `home` - customization of the home page
        - `hero.md` - customization of the front page, image from `static/img/`
        - `outline.md` - general outline of the front page, image from `static/img/`

    - `reading`
        - `_index.Rmd` - general reading assignment description
        - `01-reading.Rmd` - individual reading assignments, with links

    - `reference`
        - `_index.Rmd` - general citations and bibliography description, with BibTex file download from `static/bib`
        - `*.Rmd` - other files

    - `schedule`
        - `_index.Rmd` - schedule page, the table is rendered from `data/schedule.yaml`

    - `syllabus`
        - `_index.Rmd` - syllabus page

- `static`
    - `assignment` - working folder
    - `bib` - BibTex `references.bib` and `*.csl` files
    - `class` - working folder
    - `data` - `*.csv` files and other data
    - `files` - misc files
    - `img` - images for general pages
    - `reference` - working folder
    - `slides`
        - `image` - images for slides
        - `01-introduction_slides.Rmd` - Xarigan presentations