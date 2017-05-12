library(shiny)
library(shinyjs)
library(ggplot2)
library(ggrepel)
library(svglite)
library(colourpicker)

# Contact info footer
contact <- div(class = "text-center",
               hr(),
               "This app was made by Brian Dawes for ",
               a(href = "http://www.evanrosenlab.net/", "Evan Rosen's lab"),
               "at BIDMC.",
               br(),
               "Please send any questions or comments to ",
               a(href = "mailto:bdawes@broadinstitute.org", "bdawes@broadinstitute.org")
)

# Define app tab
app <- tabPanel("App", 
  useShinyjs(),
                
  # Some CSS styling
  tags$head(
    tags$style(HTML("
      .x-large-text {
        font-size: x-large;
      }
      .large-text {
        font-size: large;
      }"
    ))
  ),
  
  # Main controls panel
  tabsetPanel(
    tabPanel("Input",
      div(class = "well",
        # Upload, download, and plot size labels
        fluidRow(class = "x-large-text",
                 column(4, "Upload data")
        ),
        # Upload and download buttons and plot size selector
        fluidRow(
          column(4, fileInput('file1', label = NULL, accept = c(".txt", ".csv", ".tsv")))
        ),
        
        # Plot type and transparency
        fluidRow(
          class = "x-large-text",
          column(4, span("Plot type"))
        ),
        fluidRow(
          column(4, selectInput("plot_type", label = NULL,
                                choices = list("Volcano plot" = 0, "MA plot" = 1, "Abundance plot" = 2)))
        ),
        
        # Filters
        span(class = "x-large-text", "Filter data"),
        fluidRow(
          column(4, numericInput("FDR", 0.25, min = 0, max = 1, step = 0.05, label = "FDR")),
          column(4, numericInput("FC", 0, min = 0, step = 0.1, label = "logFC")),
          column(4, numericInput("CPM", 0, min = 0, step = 0.1, label = "logCPM"))
        ),
        
        # Highlight genesets
        span(class = "x-large-text", "Highlight genesets"),
        fluidRow(id = "geneset_container1",
                 column(8, textInput("geneset1", label = "Geneset", width = "100%")),
                 column(4, colourInput("color1", label = "Color", value = "#0000FF", showColour = "background", palette = "limited"))
        ),
        fluidRow(id = "change_geneset",
                 column(4, 
                        actionButton("add_geneset", "Add geneset"),
                        disabled(actionButton("rm_geneset", "Remove geneset"))
                 ),
                 column(4, radioButtons("search_style", NULL, choices = list("Exact match" = "exact", "Regex match" = "regex"), inline = TRUE))
        )
      )
    ),
    tabPanel("Output options",
      div(class = "well",
        fluidRow(class = "x-large-text",
          column(4, "Download plot"),
          column(4, "Downloaded table")
        ),
        fluidRow(
          column(4,
            div(class = "btn-group", role = "group",
              disabled(downloadButton("dl_svg", "SVG")),
              disabled(downloadButton("dl_pdf", "PDF")),
              disabled(downloadButton("dl_png", "PNG"))
            )
          ),
          column(4,
            div(class = "btn-group", role = "group",
              disabled(downloadButton("dl_csv", "CSV")),
              disabled(downloadButton("dl_tsv", "TSV"))
            )
          )
        ),
        br(),
        span(class = "x-large-text", "Plot size"),
        fluidRow(
          column(4, numericInput("width", "Width (in.)", value = 8, min = 0, step = .5)),
          column(4, numericInput("height", "Height (in.)", value = 8, min = 0, step = .5)),
          column(4, numericInput("dpi", "Resolution (dpi)", value = 300, min = 0, step = 25))
        ),
        span(class = "x-large-text", "Point options"),
        fluidRow(
          column(4, numericInput("point_size", "Size", min = 0, step = 0.5, value = 1)),
          column(4, numericInput("highlight_point_size", "Highlight Size", min = 0, step = 0.5, value = 5)),
          column(4, numericInput("alpha", "Transparency", min = 0, max = 1, step = 0.1, value = 0.3))
        ),
        fluidRow(class = "x-large-text",
          column(4, "Text size"),
          column(4, offset = 4, "Grid")
        ),
        fluidRow(
          column(4, numericInput("plot_text_size", "Plot text size", min = 0, step = 2, value = 16)),
          column(4, numericInput("gene_text_size", "Highlight text size", min = 0, step = 2, value = 5)),
          column(4, radioButtons("grids", NULL, list("Show" = "show", "Hide" = "hide"), selected = "show"))
        )
      )
    )
  ),
  # Run button
  fluidRow(
    column(4, offset = 4,
           disabled(actionButton("run", "Update plot", class="btn btn-primary btn-lg btn-block")))
  ),
  br(),
  
  # Number of up and down regulated genes/peaks
  uiOutput("up_down"),
  
  # Plot
  imageOutput("plot", height = "auto"),
  
  # Data table
  dataTableOutput("data_table"),
  
  # Add contact info
  contact
)

# Define about tab
about <- tabPanel("About",
  HTML("
<div class=\"row\">
	<nav class=\"col-md-3\">
		<ul class=\"nav nav-pills nav-stacked affix\">
			<li><a href=\"#introduction\">Introduction</a></li>
			<li><a href=\"#upload\">Upload format</a></li>
			<li><a href=\"#plotting\">Plotting</a></li>
			<li><a href=\"#filtering\">Filtering</a></li>
			<li><a href=\"#genesets\">Adding genesets</a></li>
			<li><a href=\"#table\">Table</a></li>
			<li><a href=\"#about\">About the app</a></h2>

		</ul>
	</nav>
	<div class=\"col-md-9\" style=\"margin-left: 150px\">
		<h2><a name=\"introduction\">Introduction</a></h2>
		<p>Volcano Plotter is an app that allows you to quickly and interactively generate plots from the results of differential expression analyses. Volcano Plotter can generate volcano plots, MA plots, and gene abundance plots. Significance thresholds can be set based on false discovery rate (FDR), gene abundance in log transformed counts per million (logCPM) or log transformed fold change (logFC). Important genes and sets of genes can be highlighted in the plots. Significant genes are also displayed in an interactive table. Both plots and tables can be downloaded in several common formats.</p>
		<h2><a name=\"upload\">Upload format</a></h2>
		<p>The uploaded file must be a table that can be read by R&#8217;s <a href=https://stat.ethz.ch/R-manual/R-devel/library/utils/html/read.table.html><code>read.table()</code></a> function with the options <code>header = TRUE, stringsAsFactors = FALSE</code>. The following fields must be present in the file (case sensitive):</p>
		<ul>
			<li>genes</li>
			<li>logFC</li>
			<li>logCPM</li>
			<li>PValue</li>
			<li>FDR</li>
		</ul>
		<p>Additionally, if fields beginning with Cond1_ and Cond2_ are present, these will be interpreted as logCPM values for the two conditions being compared. For example, if you have conditions WT and KO you should have fields named Cond1_WT and Cond2_KO. If these fields are not present, you will not be able to generate gene abundance plots.</p>
		<h2><a name=\"plotting\">Plotting</a></h2>
		<p>Once the data is uploaded, hit the update plot button to generate a plot. The specific type of plot is determined by the option selected from the plot type dropdown menu. The possibilities are:</p>
		<ul>
			<li><strong>Volcano plot:</strong> This plot has logFC values on the x-axis and -log10(PValues) on the y axis.</li>
			<li><strong>MA plot:</strong> This plot has logCPM values on the x-axis and logFC values on the y-axis.</li>
			<li><strong>Gene abundance plot:</strong> This plot has logCPM values for your two experimental conditions on the x- and y-axes. This option is only present if this information was included in your uploaded file (see <a href=\"#upload\">upload format</a>).</li>
		</ul>
		<p>In each plot, significant genes are highlighted in red while non-significant genes are colored gray (see <a href=\"#filtering\">filtering</a>). Additionally, you can add your own genes of interest which will be highlighted and labeled on the plot (see <a href=\"#genesets\">adding genesets</a>).</p>
		<p>Plots can be downloaded in SVG, PDF, or PNG format using the download data buttons. The downloaded file will be named according to the plot type. 
		<h2><a name=\"filtering\">Filtering</a></h2>
		<p>Volcano Plotter allows you to interactively change significance thresholds. The filters are listed under the filter data heading. The filters are:</p>
		<ul>
			<li><strong>FDR:</strong> Significant genes must be &ge; this number</li>
			<li><strong>logFC:</strong> Significant genes must be &ge; the absolute value of this number</li>
			<li><strong>logCPM:</strong> Significant genes must be &ge; this number</li>
		</ul>
		<p>Filters will not apply until the update plot button is clicked. Significant genes are highlighted in red in plots and are displayed in the table. The total number of significant genes, as well as the number of up and down regulated genes are displayed above the plot.</p>
		<h2><a name=\"genesets\">Adding genesets</a></h2>
		<p>Volcano Plotter allows you to input a gene or a list of multiple genes to highlight in the plots. Gene names are determined by the genes field of the uploaded file (see <a href=\"#upload\"> upload format</a>). Highlighted genes are colored and labeled in the plot. To enter a geneset, use the geneset field under the highlight genesets header. To enter multiple genes, separate them with commas. Spaces around the comma are ignored. For example, a geneset might look like \"IL6, IL10, IL18\".</p>
		<p>The search mode for genes can be switched between exact match and regex match. In exact match mode, the app will search for a case-insensitive exact match. In regex match mode, the app will search using case insensitive perl style regular expressions. For example, in regex match mode \"^IL\\d\\d*$\" will return every gene starting with IL followed by one or more digit and nothing else.</p>
		<p>Multiple different genesets can be entered and colored differently using the add and remove geneset buttons. The color of each geneset can be changed by clicking on the color field. By default, the first four genesets are blue, orange, green, then purple. All further genesets default to black. Any gene that is found in multiple genesets will only</p>
		<h2><a name=\"table\">Table</a></h2>
		<p>When the update plot button is clicked, an interactive table containing all the significant genes is generated below the plot. This table will contain all the fields in the input file, not just those required to run the app (see <a href=\"#upload\">upload format</a>). All the fields can be searched in the table using the top search box, while find can be searched individually using the bottom search boxes. The table can be filtered on each field by clicking on the table headers.</p>
		<p>The table can be downloaded in CSV or TSV format using the download data buttons. All significant genes will be downloaded, not just those that are currently displayed from any filters in the table.</p>
		<h2><a name=\"about\">About the app</a></h2>
		<p>Volcano plotter was made using the R shiny library. The shinyjs library is also used for some additional functionality. The colourpicker library is used to generate the color input boxes. Plotting is done using the ggplot2 library and the ggrepel library is used to position labels on highlighted genesets. Conversion to svg is done using the svglite library.</p> 
	</div>
</div>
  "),
  # Add contact info
  contact
)

# Define UI as a navbar with app and about tabs
shinyUI(navbarPage("Volcano plotter", app, about))
