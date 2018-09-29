library(shiny)
library(shinyjs)
jscode <- "shinyjs.refresh = function() { history.go(0); }"

shinyUI(fluidPage(
  tags$head(
      tags$style(HTML("
                /* The alert message box */
              .alert {
                 padding: 20px;
                 background-color: #6699ff; /* blue */
                 color: white;
                 margin-bottom: 15px;
              }
              .alert_red {
                 padding: 20px;
                 background-color: #f44336; /* red */
                 color: white;
                 margin-bottom: 15px;
              }
              .alert_warning {
                 padding: 20px;
                 background-color: #00cc44; /* green */
                 color: white;
                 margin-bottom: 15px;
              }
                 "))
      ),

  titlePanel("The effectR package graphical user interface."),
    mainPanel(
      tags$hr(),
      h4("effectR predicts effectors for genomes in a fast and reproducible way."),
      tags$hr(),
      singleton(tags$body(HTML(
        '
        <div class="alert">
        <h3> Instructions and guidelines:</h3>
        <b> Please read the instructions very carefully. Follow each step and be patient, the processes might take some time to start. CLICK ONLY ONE IN EACH BUTTON.</b>
        <hr>
        <p></p>
To obtain a better assesment of the effectors in your genome please use a <b>6-frame translation</b> of the genome of interest (Use <b><a href="http://emboss.sourceforge.net/apps/cvs/emboss/apps/getorf.html">getorf</a></b> to obtain the 6-frame translation and use the output file as input in this pipeline.)
<hr>
<b>Example file:</b> Subset of Twenty-five sequences of <i>Phytophthora infestans</i> obtained from <a href=http://fungidb.org/>FungiDB</b></a> and processed in <b>getorf</b>
<p></p>
<p></p>
<a href="https://www.dropbox.com/s/td4edffnxsdmu7s/test_infestans.fasta?dl=1" class="btn btn-warning btn-med" role="button"><span class="glyphicon glyphicon-floppy-saved" aria-hidden="true"></span> Example dataset</a>
        </div>
        '
        ))),
      fileInput('file1', 'Upload FASTA formatted files (.fsa,.fasta,.fas,.orf).',accept=c('text','.fsa','.fasta','.fas','.orf')),
      verbatimTextOutput('input_fasta'),
      singleton(tags$body(HTML(
        '
        <div class="alert_red">
        <b>WARNING:</b> Please wait until the prompt returns the number of sequences processed to proceed to the REGEX steps. The process might take minutes, <b>please be patient</b>.
        </div>
        '
      ))),
      useShinyjs(),
      extendShinyjs(text = jscode),
      actionButton("refresh", "Restart app", style="color: #fff; background-color: #b72103; border-color: #a31f04"),


      tags$hr(),
      h2("Step 1"),
      h3("Regular expression search"),
      p("Step 1 uses a regular expression search to identify the sequences with the motifs of interest."),
      p(),
      selectInput("motif_sel", label = h3("Select motif to search"),
                  choices = list("RxLR" = "RxLR", "CRN" = "CRN"), selected = "RxLR"),
      p(),
      actionButton('regex_seach', "Do REGEX search", icon("paper-plane"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
      p(),
      textOutput('regex_wait'),
      verbatimTextOutput('regex'),
      singleton(tags$head(HTML(
        '
        <script type="text/javascript">
        $(document).ready(function() {
        // disable download at startup. downloadHMM is the id of the downloadButton
        $("#downloadREGEX").attr("disabled", "true").attr("onclick", "return false;");

        Shiny.addCustomMessageHandler("regex_ready", function(message) {
        $("#downloadREGEX").removeAttr("disabled").removeAttr("onclick").html(
        "<i class=\\"fa fa-download\\"></i>Download (file size: " + message.fileSize + ")");
        });
        })
        </script>
        '))),
      downloadButton('downloadREGEX', 'Download REGEX candidates'),
      helpText("Download will be available once the processing is completed."),

      #ADD A DOWNLOADHANDLER TO OBTAIN THE REGEX SEQUENCES

      tags$hr(),
      h2("Step 2"),
      h3("HMM search"),
      p("Step 2 aligns the candidate genes obtained in Step 1 to create a Hmmer search to identify additional potential effectors."),
      p(),
      textInput("seed", "Random seed number. (This will affect the number of predicted proteins by HMMER. Use the same seed for reproducible results).",value = '12345'),
      actionButton("do", "Align and do HMM search", icon("paper-plane"), style="color: #fff; background-color: #61c66d; border-color: #21a831"),
      p(),
      verbatimTextOutput('mafft'),
      verbatimTextOutput('hmmer_press'),
      verbatimTextOutput('hmmer_search'),
      # JAVASCRIPT: Download button
      singleton(tags$head(HTML(
        '
        <script type="text/javascript">
        $(document).ready(function() {
        // disable download at startup. downloadHMM is the id of the downloadButton
        $("#downloadHMM").attr("disabled", "true").attr("onclick", "return false;");

        Shiny.addCustomMessageHandler("download_ready", function(message) {
        $("#downloadHMM").removeAttr("disabled").removeAttr("onclick").html(
        "<i class=\\"fa fa-download\\"></i>Download (file size: " + message.fileSize + ")");
        });
        })
        </script>
        '))),
      downloadButton('downloadHMM', 'Download HMM files'),
      helpText("Download will be available once the processing is completed."),

      h3("HMM Logo"),
      p("Create a bar plot resembling a \"sequence logo\" of the sequence alignment used to build the HMM model in Step 2."),
      actionButton("logo", "Create LOGO",icon("paper-plane"), style="color: #fff; background-color: #61c66d; border-color: #21a831"),
      p(),
      plotOutput('preImage', width = "200%"),
      singleton(tags$head(HTML(
        '
        <script type="text/javascript">
        $(document).ready(function() {
        // disable download at startup. downloadHMM is the id of the downloadButton
        $("#logo_img").attr("disabled", "true").attr("onclick", "return false;");

        Shiny.addCustomMessageHandler("logo_ready", function(message) {
        $("#logo_img").removeAttr("disabled").removeAttr("onclick").html(
        "<i class=\\"fa fa-download\\"></i>Download (file size: " + message.fileSize + ")");
        });
        })
        </script>
        '))),
      #downloadButton('logo_img', 'Download HMM alignment LOGO'),

      tags$hr(),
      h2("Step 3"),
      h2("Tables and statistics"),
      h3("Summary table"),
      p("Combine the candidate effectors in Step 1 and Step 2, remove reduntant sequences and detect the positions of the motifs of interest (RxLR and EER)"),
      actionButton("combine", "Combine datasets and return table", icon("paper-plane"), style="color: #fff; background-color: #e55e2d; border-color: #d8440f"),
      p(),
      tableOutput('final_table'),
      p(),
      h3("Motif table"),
      p("Returns a table based on Haas et al., (2009) with aminoacid sequence, positions of the MOTIFs of interest and a summary of presence/absence of motifs of interest"),
      actionButton("motif", "Motif Table", icon("paper-plane"), style="color: #fff; background-color: #e55e2d; border-color: #d8440f"),
      p(),
      tableOutput('motif_table'),
      # JAVASCRIPT: Download button
      singleton(tags$head(HTML(
        '
        <script type="text/javascript">
        $(document).ready(function() {
        // disable download at startup. downloadHMM is the id of the downloadButton
        $("#motifDownload").attr("disabled", "true").attr("onclick", "return false;");

        Shiny.addCustomMessageHandler("motif_ready", function(message) {
        $("#motifDownload").removeAttr("disabled").removeAttr("onclick").html(
        "<i class=\\"fa fa-download\\"></i>Download (file size: " + message.fileSize + ")");
        });
        })
        </script>
        '))),
      downloadButton('motifDownload', 'Download MOTIF table'),
      p(),
      h3("Motif summary table"),
      p("Table with a summary of the number of sequences with all, one or no motifs."),
      actionButton("motif_summ", "Summary Table", icon("paper-plane"), style="color: #fff; background-color: #e55e2d; border-color: #d8440f"),
      tableOutput('motif_summ_table'),
      p(),
      tags$hr(),
      h2("Download final FASTA candidate file"),
      h4("This FASTA file contains all sequences found in both steps."),
      singleton(tags$body(HTML(
        '
        <div class="alert_warning">
        <b>DISCLAIMER:</b>
        <p></p>
       We recommend a final <b>manual curation</b> and <b>verification</b> of the effectors obtained. Other non-canonical RxLR effectors might not be called due to the REGEX + HMM based searches.
        <p></p>
        In addition, we recommend the use of <b><a href="http://www.cbs.dtu.dk/services/SignalP-2.0/">SignalP 2.0</a></b> in the resulting output to determine if the called effectors have <b>signal peptide evidence</b>
    </div>
        '
      ))),
      downloadButton('downloadData', 'Download'),
      useShinyjs(),
      extendShinyjs(text = jscode)
    ),
    tags$hr()
)
)

