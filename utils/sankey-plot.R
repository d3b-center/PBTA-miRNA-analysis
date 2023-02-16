# Sankey Plot
# Author: Shehbeel Arif
# Children's Hospital of Philadelphia

# Load library
library(plotly)

# Build Sankey plot
fig <- plot_ly(
  # Specify Sankey plot
  type = "sankey",
  # Display plot horizontally
  orientation = "h",
  
  # Create nodes
  node = list(
    label = c("PBTA", "ATRT", "Cranio", "Ependy", "Gangli", "GNT", "HGAT", "LGAT", "Medullo", "Schwan", "Tert",
              "G1", "G2"),
    color = c("black", "blue", "red", "green", "purple", "orange", "lightblue", "pink", "lightgreen", "violet", "yellow"),
    pad = 15,
    thickness = 20,
    line = list(
      color = "black",
      width = 0.5
    )
  ),
  # Specify links between nodes
  link = list(
    # source nodes
    source = c(0,0,0,0,0,0,0,0,0,0,4,4),
    # target nodes
    target = c(1,2,3,4,5,6,7,8,9,10,11,12,13), 
    # size of links
    value =  c(17,21,36,19,3,34,102,35,2,1,10,9) # size of link
  )
)

# Specify plot layout
fig <- fig %>% layout(
  # Plot title
  title = "Pediatric Brain Tumor Subtyping by miRNA",
  # Font size
  font = list(
    size = 10
  )
)

# Display figure
fig
