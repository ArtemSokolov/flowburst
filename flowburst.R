## Visualization of Greg's cell type signatures
##
## by Artem Sokolov

library( ggsunburst )
library( tidyverse )

## Main driver of the code in this file
main <- function()
{
    fn <- "alpha_1.0%_vectors.csv"

    ## ########## ##
    ## Pre-define ##
    ## ########## ##
    
    ## Markers
    vMarkers <- c( "B220", "CD11b", "CD11c", "CD3e", "CD4",
                  "CD45", "CD49b", "CD8a", "F480", "Ly6C", "Ly6G" )

    ## Canonical cell states
    vCanonical <- c( "B", "CD8aposT", "ISPT", "CD4posT", "NK", "Non-immune",
                    "Precursor", "Mac", "Mono", "DPT", "DNT", "PMN", "DC" )

    ## Mapping to short labels
    mClean <- c( "CD8aposT" = "CD8T", "CD4posT" = "CD4T", "CD8aposISPT" = "ISPT" )

    ## Order on the figure (to avoid label overlap)
    ## Must use "clean"/short labels from mClean
    vOrder <- c( "B", "CD4T", "DC", "CD8T", "ISPT", "Precursor", "DNT",
                "NK", "Mono", "Mac", "PMN", "DPT" )

    ## ######## ##
    ## Pipeline ##
    ## ######## ##
    
    ## Load cell type definitions and annotate each cell type with the corresponding
    ##  canonical label
    CT <- loadCellTypes( fn, vMarkers, vCanonical, mClean )

    ## Identify the set of markers to display next to each cell type
    MK <- getMarkers( CT )
    
    ## Prepare node data for the sunburst plot
    SB <- prepareNodeData( CT, vOrder ) %>% genMarkerLbls( MK )

    ## Plot the file result
    plotFlowburst( SB, "cell_types.pdf" )
}

## Given a vector and a partial map, creates a full mapping
## v - vector of values (duplicates are allowed and removed internally)
## pm - partial map of the form c( "Old Value1" = "New Value1", "OV2" = "NV2", ...)
##   "Old Value"s may or may not occur in v
## do.map - FALSE -> returns map; TRUE -> returns result of map applied to v
##   Default: TRUE
partialMap <- function( v, pm, do.map=TRUE )
{
    m <- unique(v) %>% set_names
    vCommon <- intersect( m, names(pm) )
    if( length(vCommon) > 0 ) m[vCommon] <- pm[vCommon]
    if( do.map ) unname(m[v]) else m
}

## Loads the cell subtype signatures and assigns cell type labels to each
## fn - file name of the input data
## vMarkers - vector of pre-defined flow markers
## vCanonical - vector of canonical state names; these are tail-matched against
##   individual entries in cell_type column from the end of the string
## mClean - (optional) mapping of names from vCanonical and cell_type to
##   clean/short versions that are displayed on the plot
loadCellTypes <- function( fn, vMarkers, vCanonical, mClean = NULL )
{
    X <- read_csv( fn, col_types = cols() ) %>%
        select( cell_type, one_of(vMarkers) )

    ## Match up cell types to their canonical state names
    X$Canonical <- set_names( vCanonical ) %>% map( ~grep(str_c(.x,"$"), X$cell_type) ) %>%
        reshape2::melt() %>% arrange( value ) %>% .$L1

    ## Map canonical names to their "clean" versions (as applicable)
    ## Reorder columns to show associated canonical state names first
    ## Remove non-immune from consideration because all markers are 0
    X %>% mutate( cell_type = partialMap( cell_type, mClean ),
                 Canonical = partialMap( Canonical, mClean ) ) %>%
        select( Canonical, cell_type, everything() ) %>%
        filter( Canonical != "Non-immune" )
}

## Appends a suffix to each element of v
## Returns empty vector if v is empty
str_cc <- function( v, sfx )
{
    if( length(v) == 0 ) return(character(0))
    str_c( v, sfx )
}

## Annotate each Cell Type with its markers
## Distinguish canonical cell types, for which the entire set is displayed,
##   from other cell types, for which only the difference in markers is shown
getMarkers <- function( CT )
{
    CT %>% filter( Canonical != "Non-immune" ) %>%
        ## Retrieve the list of markers for each cell type
        group_by( Canonical, cell_type ) %>%
        do( ctm = colnames(.)[.==1] ) %>% ungroup() %>%
        ## Annotate with the corresponding canonical markers
        group_by( Canonical ) %>%
        mutate( cnm = ctm[ Canonical == cell_type ] ) %>% ungroup() %>%
        ## Find the set of markers that differs between cell type and canonical
        mutate( mpos = map2(ctm,cnm,setdiff ), mneg = map2(cnm,ctm,setdiff) ) %>%
        ## Annotate markers with ^{+} and ^{-} respectively and combine into a single vector
        mutate( mpos = map( mpos, str_cc, "^'+'" ), mneg = map( mneg, str_cc, "^'-'" ) ) %>%
        mutate( mdiff = map2( mpos, mneg, c ) ) %>%
        ## Isolate the markers to be plotted
        ## This is column cnm (==ctm) for canonical types and column mdiff for all others
        mutate( markers = ifelse( Canonical == cell_type, cnm, mdiff ) ) %>%
        select( Canonical, cell_type, markers ) %>% unnest( markers )
}

## Sets up the initial sunburst structure
prepareNodeData <- function( CT, vOrder )
{
    ## Annotate internal nodes with corresponding canonical cell types, which is
    ##   basically themselves. This is for coloring purposes.
    XCn <- data_frame( node = unique(CT$Canonical) ) %>% mutate( Canonical = node ) %>%
        slice( match(vOrder, node) )

    ## Generate node-parent associations for non-canonical cell types
    XX <- CT %>% select( node = cell_type, Canonical ) %>% mutate( parent = Canonical ) %>%
        filter( !(node == parent) ) %>% bind_rows( XCn )

    ## Store the result to a temp file to be read by the Python script
    tf <- tempfile()
    XX %>% write_csv( tf, na="" )
    sb <- sunburst_data( tf, type = "node_parent", sep=",", node_attributes = "Canonical" )
    
    ## Adjust rectangle height
    sb$rects <- sb$rects %>% mutate( ymin = ifelse( ymin == min(ymin), ymin + 0.5, ymin ),
                                    ymax = ifelse( ymax == max(ymax), ymax - 0.75, ymax ) )
    sb$node_labels <- sb$node_labels %>% mutate( y = y + 0.25 )
    sb$leaf_labels <- sb$leaf_labels %>% mutate( y = ifelse( y_out < 0, y + 0.25, y - 0.625 ) ) %>%
        mutate( yend = ifelse( Canonical == "DC", max(y)+.325, max(y)+.6 ) )
    sb
}

## Generates marker annotations for inner and outer rings
## Inputs:
##   SB - sunburst object, as created by prepareNodeData
##   MK - marker associations, as computed by getMarkers
##
## Returns modified SB object with two new fields:
##   LL - outer/leaf marker labels
##   NN - inner/node marker labels
genMarkerLbls <- function( SB, MK )
{
    ## Create external/leaf annotations
    SB$LL <- SB$leaf_labels %>% rename( cell_type = label ) %>%
        mutate_at( vars(cell_type, Canonical), as.character ) %>%
        inner_join( MK ) %>%
        mutate( y = ifelse( cell_type == "DC", max(y) + .425, max(y) + .7 ) ) %>%
        group_by( cell_type ) %>% mutate( y = y + (1:length(y) - 1) * 0.125 ) %>% ungroup
    
    ## Create node annotations
    SB$NN <- SB$node_labels %>% rename( cell_type = label ) %>%
        mutate_at( vars(cell_type, Canonical), as.character ) %>%
        inner_join( MK ) %>%
        mutate( y = min( SB$rects$ymin ) - 0.1 ) %>% group_by( cell_type ) %>%
        mutate( y = y - (1:length(y)-1) * 0.125 ) %>% ungroup

    SB
}

## Plots the final sunburst object, after it's been passed through genMarkerLbls()
plotFlowburst <- function( SB, fnOut = "cell_types.pdf" )
{
    ## Define the palette
    v <- SB$rects %>% select( name, Canonical ) %>% mutate_all( as.character ) %>%
        filter( name == Canonical ) %>% magrittr::extract2( "name" )
    pal <- scales::hue_pal()(length(v))[match( sort(v), v )]

    ## Plot everything
    gg <- sunburst( SB, rects.fill.aes = "Canonical",
                   node_labels.min = 1, node_labels.size = 5,
                   leaf_labels = FALSE, node_labels = FALSE ) + guides( fill = FALSE ) +
        scale_fill_manual( values = pal ) +
        geom_segment( data = SB$leaf_labels, aes( x = x, xend=x, y=y+.375, yend=yend ) ) +
        geom_text( aes( x=x, y=y, angle=ifelse( label == "CD4T", rangle, pangle ), label=label ),
                  size = 5, data = SB$node_labels, fontface="bold" ) +
        geom_text( aes( x=x, y=y, angle=angle, label=label, hjust=0.5 ), size=4,
                  data = subset(SB$leaf_labels, y_out != 0), fontface="bold" ) +
        geom_text( aes( x=x, y=y, label=markers, angle=pangle, hjust=0.5 ),
                  size=4.5, data=SB$LL, parse = TRUE ) +
        geom_text( aes( x=x, y=y, label=markers, angle=pangle, hjust=0.5 ),
                  size=4.5, data=SB$NN ) +
        scale_y_continuous( limits=c(-2.9, -0 ), expand=c(0,-0.55) ) +
        theme( plot.margin = unit( c(0,0,0,0), "mm" ),
              axis.title.x = element_blank(), axis.title.y = element_blank() )
    
    ggsave( fnOut, gg, width=8, height=8 )
}

