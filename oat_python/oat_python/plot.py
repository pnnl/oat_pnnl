#   Visualizatino of simplicial complexes
#
#   Greg Henselman-Petrusek Roek
#   2022-09-27

import oat_python
import oat_python.barcode

import itertools
import plotly
import plotly.express

import plotly.graph_objects as go
import numpy as np
from . import hypergraph
from cmath import inf





#   ========================================
#   PERSISTENCE DIAGRAM + BARCODE
#   ========================================





def bdtpairs__trace( birth_death_time_pairs , upper_limit = None ):
    """
    Plot the barcode.  Each bar of form (birth, infinity) will be
    represented by a point of form (birth, upper_limit).

    :param birth_death_pairs: an iterable of (birth,death) pairs
    :param upper_limit: a real number
    :returns: a plotly trace of type go.Scatter
    """

    x = [ x[0] for x in birth_death_time_pairs ]
    y = [ np.minimum( x[1], upper_limit ) for x in birth_death_time_pairs ]
    trace = go.Scatter( x=x, y=y, mode='markers' )
    return trace

def pd_guidelines( upper_limit ):
    """
    Given a scalar value `u = upper_limit`, returns
    - a trace for the line segment (0,u) -- (u,u); dashed
    - a trace for the line segment (0,0) -- (u,u); solid
    """
    traceh = go.Scatter(x=[0,upper_limit], y=[upper_limit,upper_limit], mode="lines")
    traceh.update( line=dict( dash="dot", color="black" ) )
    traced = go.Scatter(x=[0,upper_limit],y=[0,upper_limit], mode="lines", )
    traced.update( line=dict( color="black" ) )    
    return traced, traceh

def pd( df=None, guideline_limit = None ):
    """
    Return figure for the persistence diagram.
    :param df: a dataframe with columns
    :param guidline_limit: sets the limit for the diagonal (=1.2 * guideline_limt) and horizontal (=1.1 * guideline_limit)
    visual guidelines in the diagram

    The hovertext in the resulting diagram provides a "Feature id number" for each point `P`.  This is a unique id that
    can be used to 
    """
    fig                         =   go.Figure( data=[] )


    from cmath import inf  
    C                               =   guideline_limit
    if C is None:
        finite_endpoints = { x for x in df["birth"].tolist() + df["death"].tolist()  if x != inf }
        if len(finite_endpoints)==0:
            C                       =   1
            infinity_proxy          =   1 # determines whwere we will draw the "infinity" line
            diagonal_limit          =   1.1 # how much farther than "infinity proxy" 
        else:
            fin_max                 =   max( finite_endpoints )
            fin_min                 =   min( finite_endpoints )
                        
            C                       =   max( finite_endpoints )
            if fin_max > fin_min:
                infinity_proxy      =   C + 0.1 * (fin_max - fin_min)
                diagonal_limit      =   C + 0.2 * (fin_max - fin_min)
            else:
                infinity_proxy      =   C + 1.1 * np.abs(C)
                diagonal_limit      =   C + 0.2 * np.abs(C)

    


    fig.add_hline(y= infinity_proxy, line=dict(dash="dot"))    
    trace_lined, trace_lineh    =   pd_guidelines( diagonal_limit ) 
    trace_lined.update(showlegend=False)#name="x=y")
    fig.add_trace(trace_lined)    

    import plotly.express as px

    # fig = px.scatter()
    colors = px.colors.qualitative.Plotly;

    
    flagged_dimensions  =   set()

    for index, row in df.iterrows():
        x               =   row["birth"]
        y               =   np.minimum( row["death"], infinity_proxy )
        dimension       =   row["dimension"]
        color           =   colors[ dimension ]
        text            =   f"birth filtration = {x}<br>" + \
                            f"death filtration = {y}<br>" +\
                            f"interval length = { row['death'] - row['birth']}<br>" +\
                            f"birth simplex = {row['birth simplex']}<br>"
        
        if 'cycle nnz' in row.keys():
            text += f"cyle representative nnz = {row['cycle nnz']}<br>"
        if 'bounding nnz' in row.keys():
            text += f"bounding chain nnz = {row['bounding nnz']}<br>"

        text            +=  f"row of dataframe = {index}"
        
        trace   =   go.Scatter(
                        x                   =   [x],
                        y                   =   [y],
                        mode                =   "markers",
                        text                =   [text],
                        showlegend          =   dimension not in flagged_dimensions,
                        legendgroup         =   dimension,
                        marker              =   dict(color=color),
                        name                =   f"Dimension {dimension}",
                    )
        fig.add_trace(trace)

        flagged_dimensions.add(dimension)
        

    # for value in df['dimension'].unique():
    #     color = colors[ value % len(colors) ]


    #     dimslice = df[df['dimension'] == value].copy() # makes a copy
    #     dimslice.drop("cycle representative", axis=1)
    #     dimslice.drop("bounding chain", axis=1)
    #     dimslice.replace( inf, infinity_proxy, inplace=True) # make infinite points finite)
    #     dimslice.rename( columns={"birth":"birth filtration", "death":"death filtration", "cycle nnz":"cycle representative nnz", "bounding nnz":"bounding chain nnz"}, inplace=True )
    #     dimslice["birth simplex"] = [ " " + str(x) for x in dimslice["birth simplex"]]
    #     dimslice["death simplex"] = [ " " + str(x) for x in dimslice["death simplex"]]        
    #     dimslice["row of data frame"] = dimslice.index
    #     trace = px.scatter(dimslice, x="birth filtration", y="death filtration", hover_data=["birth filtration", "death filtration", "birth simplex", "death simplex", "cycle representative nnz", "bounding chain nnz", "row of data frame"],).data[0]
    #     trace.update( name = f"Dimension {value}", showlegend=True, marker=dict(color=color))
    #     fig.add_trace( trace )

    # fig.show()    

    # dims = list( { bar.dimension() for bar in barcode.bars() } )
    # dims.sort()

    # for dim in dims:
    #     intervals               =   barcode.intervals(dim)
    #     births                  =   [x[0] for x in intervals]
    #     deaths                  =   [np.minimum(x[1], 1.1*C) for x in intervals]
    #     hovertext               =   [f"birth {x}<br>death {y}<br>homology dimension {dim}<br>feature id number {id}" \
    #                                  for (x,y,id) in intervals ]
    #     trace                   =   go.Scatter( x=births, y=deaths, mode='markers' )
    #     trace.update(hoverinfo="text")
    #     trace.update(name=f"Dim {dim}")
    #     trace.update(hovertext=hovertext)
    #     fig.add_trace(trace)

    fig.update_layout(yaxis_range=[-0.1 * C, 1.2 * C ].sort())    
    fig.update_layout(xaxis_range=[-0.1 * C, 1.2 * C ].sort())    
    fig.update_layout(height=600,width=650)    
    # fig.update_layout(scene=dict(aspectmode='manual', aspectration=go.layout.Aspectraio()))
    return fig



def barcode( barcode=None, guideline_limit = None ):
    """
    Returns a Plotly figure for the barcode of a filtered chain complex.

    :param barcode: an Pandas data frame with `birth`, `death`, and `dimension` columns.
    :param guidline_limit: sets the limit for the vertical (=1.1 * guideline_limt) visual guidelines in the diagram

    The hovertext in the resulting diagram provides a "Feature number" for each point `P`.  In the feature number is `m`
    and the Feature dimension is `d`, then the `m`th element in the list of dimension `d` cycle representatives is the
    cycle representative for point `P`.
    """
    fig                         =   go.Figure( data=[] )

    # fig.add_trace( go.Scatter(x=[],y=[],yaxis="y2") )   
    C                           =   guideline_limit
    if C is None:
        C                       =   oat_python.barcode.max_finite_value( barcode["birth"].tolist() + barcode["death"].tolist() )
        if C is None:
            C                   =   1

    fig.add_vline(x= 1.1* C, line=dict(dash="dot"))    


    intervals                   =   sorted(
                                        list(
                                            zip(
                                                barcode["dimension"],
                                                barcode["birth"], 
                                                barcode["death"], 
                                                barcode.index 
                                            )
                                        )
                                    )    
    
    color_sequence              =   plotly.express.colors.qualitative.Plotly
    num_colors                  =   len(color_sequence)
    max_dim                     =   max(barcode["dimension"])
    x                           =   [[] for _ in range(max_dim + 1)]
    y                           =   [[] for _ in range(max_dim + 1)]
    hovertext                   =   [[] for _ in range(max_dim + 1)]

    # coordinates and hover text for bars
    for bar_counter, (dim,birth,death,id) in enumerate(intervals):
        x[dim].append( birth )
        x[dim].append( np.minimum( death, 1.1 * C ) )
        x[dim].append( np.nan )
        y[dim].append( bar_counter )
        y[dim].append( bar_counter )
        y[dim].append( np.nan )    
        newtext                 = \
        "birth filtration {birth}" \
        + f"death filtration {death}" \
        + f"birth simplex {barcode['birth simplex'][id]}" \
        + f"death simplex  {barcode['death simplex'][id]}" \
        + f"cycle representative nnz {barcode['cycle nnz'][id]}" \
        + f"bounding chain nnz {barcode['bounding nnz'][id]}" \
        + f"row of data frame {id}"
        hovertext[dim].append( f"birth {birth}<br>death {death}<br>homology dimension {dim}<br>row of dataframe {id}" )         

    for dimension, (x,y,hovertext) in enumerate(zip(x,y,hovertext)):
        # traces for bars
        color                   =   color_sequence[ dimension % num_colors ]        
        trace                   =   go.Scatter( x=x, y=y, name=f"Dim {dimension}", mode='lines', line=dict(color=color), )
        hov                     =   [ t for t in hovertext for _ in range(3) ]
        trace.update(hoverinfo="text")
        trace.update(hovertext=hov)
        fig.add_trace(trace) 

        # # traces for hover text
        # trace                   =   go.Scatter( x=x[::3], y=y[::3], mode='markers+text', marker=dict(color=color), legendgroup=dimension, )        
        # trace.update(hoverinfo="text")
        # trace.update(name=f"Dim {dimension}")
        # trace.update(hovertext=hovertext)

        # add trace to figure
        # fig.add_trace(trace)     

    fig.update_layout(yaxis_range=[ 0, bar_counter + 1 ])    
    fig.update_layout(xaxis_range=[-0.1 * C, 1.2 * C ])    
    fig.update_layout(height=600,width=600)    
    return fig




#   =======================
#   SCATTER
#   =======================


#   DEPRECATED
#
# def ez_scatter3( coords_as_rows=np.zeros((0,3)), color = None, name="", cmin=0, cmax=1):
#     """
#     Convenience function for plotting point clouds in R^3
#     """
#     print("The developers are considering deprecating this function.")
#     coords_as_rows = np.array(coords_as_rows)
#     if color is None:
#         return go.Scatter3d( x = coords_as_rows[:,0], y=coords_as_rows[:,1], z=coords_as_rows[:,2], mode="markers", marker=dict(symbol="circle-open"), name=name, showlegend=True, )    
#     else:
#         return go.Scatter3d( x = coords_as_rows[:,0], y=coords_as_rows[:,1], z=coords_as_rows[:,2], marker=dict(symbol="circle-open", color=color, cmin=cmin, cmax=cmax, colorscale="Viridis"), mode="markers", name=name, showlegend=True, )


#   =======================
#   MDS FOR PLOTTING
#   =======================


#   HELPER FUNCTIONS TO GENERATE MDS HOP COORDINATES
#   -------------------------------------------------

from sklearn import manifold
import networkx as nx
import networkx as nx
import itertools

def adjacency_from_simplices(simplices):
    """
    simplices: a list of length-k lists
    """
    import itertools
    vertices = [ v for simplex in simplices for v in simplex  ]
    num_node = np.max(vertices) + 1 
    
    A = np.zeros((num_node, num_node))
    for simplex in simplices:
        for i,j in itertools.combinations(simplex,2):
            A[j, i] = 1
            A[i, j] = 1
            A[i, i] = 1
            A[i, j] = 1
    return A

def nx_graph_from_simplices(simplices):
    """
    Generate an nx graph representing the 1-skeleton of the simplicial complex.

    :param  simplices - a list of length-k lists
    :return: an nx graph representing the 1-skeleton of the simplicial complex
    """
    G=nx.Graph()
    
    for simplex in simplices:
        for v in simplex:
            G.add_node(v)
        for i,j in itertools.combinations(simplex, 2):
            G.add_edge(i,j)
    return G

def hop_distance_from_nx_graph(G):
    hop_distances = dict(nx.shortest_path_length(G))
    n_vertices = G.number_of_nodes()
    vertex_labels = [x for x in hop_distances.keys() ]
    D  = np.zeros( (n_vertices, n_vertices) )
    for vertex in range(n_vertices):
        for neighbor in range(vertex,n_vertices):
            D[vertex][neighbor] = hop_distances[ vertex_labels[vertex]  ][ vertex_labels[neighbor] ]
            D[neighbor,vertex] = D[vertex,neighbor]
    return D, vertex_labels


def hop_mds_from_nx_graph(G, dimension=3):
    D, vertex_labels    = hop_distance_from_nx_graph(G)

    mds = manifold.MDS(
        n_components    =   dimension,
        dissimilarity   =   "precomputed",
    )
    coords = mds.fit(D).embedding_    

    return { vertex_labels[k]: coords[k] for k in range(len(vertex_labels)) }



def hop_mds_from_simplices(simplices, dimension=3):
    """
    Returns a dictionary of length-3 1-d arrays
    
    :param simplices: a list of length-k lists
    :param dimension: the dimension of the space were we wish to embed our points
    """
    G = nx_graph_from_simplices(simplices)
    return hop_mds_from_nx_graph(G, dimension=dimension)





#   ========================================
#   EDGES
#   ========================================



def edge__trace2d( edge, coo ):
    """
    Generate a 2d trace for the given edge

    :param edge: an sequence (e.g. a list, array, etc.) of two integers, representing the vertices
    :param coo:  an indexible object that maps integers to x-y-z coordinates, e.g. a dictionary of tuples

    :return trace: a Plotly Scatter trace for the edge, with an "extra" vertex in the center (which is invisible but
    allows the user to view over data when the cursor passes over the center point)

    **For developers** It would be possible to offer a batched function with takes a sequences of
    edges and returns a single trace (separating individual edges with `None` entries).  We don't
    currently offer this because (1) returning traces edge-by-edge gives a practical advantage in terms
    of fine-grained control, in several common use cases, and (2) one can simulate the behavior of
    a single trace using the `legendgroup` keyword from Plotly.    
    """
    x0, x1  =   coo[edge[0]][0], coo[edge[1]][0]
    y0, y1  =   coo[edge[0]][1], coo[edge[1]][1]
    x       =   [ x0, (x0+x1)/2, x1 ]
    y       =   [ y0, (y0+y1)/2, y1 ]
    return go.Scatter(x=x, y=y, mode="lines")

    # x       =   []
    # y       =   []
    # data    =   []
    # for edge in edges:
    #     x = x + [ coo[edge[0]][0], coo[edge[1]][0], np.nan ]
    #     y = y + [ coo[edge[0]][1], coo[edge[1]][1], np.nan ]        
    #     data.append( go.Scatter(x=x,y=y, showlegend=False) )
    # return data


def edge__trace3d( edge, coo ):
    """
    Generate a 2d trace for the given edge

    :param edge: an sequence (e.g. a list, array, etc.) of two integers, representing the vertices
    :param coo:  an indexible object that maps integers to x-y-z coordinates, e.g. a dictionary of tuples

    :return trace: a Plotly Scatter3d trace for the edge, with an "extra" vertex in the center (which is invisible but
    allows the user to view over data when the cursor passes over the center point)

    **For developers** It would be possible to offer a batched function with takes a sequences of
    edges and returns a single trace (separating individual edges with `None` entries).  We don't
    currently offer this because (1) returning traces edge-by-edge gives a practical advantage in terms
    of fine-grained control, in several common use cases, and (2) one can simulate the behavior of
    a single trace using the `legendgroup` keyword from Plotly.
    """
    x0, x1  =   coo[edge[0]][0], coo[edge[1]][0]
    y0, y1  =   coo[edge[0]][1], coo[edge[1]][1]
    z0, z1  =   coo[edge[0]][2], coo[edge[1]][2]    
    x       =   [ x0, (x0+x1)/2, x1 ]
    y       =   [ y0, (y0+y1)/2, y1 ]
    z       =   [ z0, (z0+z1)/2, z1 ]        
    return go.Scatter3d(x=x, y=y, z=z, mode="lines")



def edges__trace3d( edges, coo ):
    """
    Generate a 2d trace for the given edges

    :param edges: an iterable that iterates over pairs of vertices (tuples, lists, etc.)
    :param coo:  an indexible object that maps integers to x-y-z coordinates, e.g. a dictionary of tuples    

    **Note for developers** This method tends to be less flexible than repeated calls to `edge__trace3d`,
    however the developers have kept it in the library thus (1) to preserve a well-tested function, and (2)
    to conserve at least one functional example of how to plot multiple disjoint edges in a single trace.
    """
    x       =   []
    y       =   []
    z       =   []
    data    =   []
    for edge in edges:
        x = x + [ coo[edge[0]][0], coo[edge[1]][0], None ]
        y = y + [ coo[edge[0]][1], coo[edge[1]][1], None ] 
        z = z + [ coo[edge[0]][2], coo[edge[1]][2], None ]                
    return go.Scatter3d(x=x,y=y,z=z)




#   ============================
#   TRIANGLES
#   ============================

def triangle__trace3d(triangle, coo=[]):
    """
    Generates a Plotly `Mesh3d` trace.

    :param triangles_as_rows is an array or list of lists of
    integers, where each integer points to a row of 
    :param coo: an indexible object that maps
    integers to x-y-z coordinates, e.g. a dictionary of tuples  
    """

    # Plotly Mesh3d requires x, y, z to be formatted as lists, tuples, or arrays
    # for this reason, to avoid bundling a large point cloud into this trace, we
    # relabel the vertices incident to the simplices with labels {1, ..., n}, and
    # then generate coordinate array with only as many coordinates as we need
    triangle = np.array(triangle)
    nvl2ovl = np.unique( np.array(triangle).flatten() )
    ovl2nvl = { nvl2ovl[nvl]: nvl for nvl in range(len(nvl2ovl))}
    triangle_relabeled = [ ovl2nvl[x] for x in triangle ]

    coordinates_as_rows = np.array( [ coo[ v ] for v in nvl2ovl ] )

    trace   =   go.Mesh3d(
                    x=coordinates_as_rows[:,0],
                    y=coordinates_as_rows[:,1],
                    z=coordinates_as_rows[:,2],
                    i = [triangle_relabeled[0]],
                    j = [triangle_relabeled[1]],
                    k = [triangle_relabeled[2]],
                )    
    return trace 


def triangle__trace2d( triangle=[], *, coo=[]):
    """
    Generates a filled plotly trace with `go.Scatter(x=x, y=y, fill="toself",)`
    :param triangle: a length-3 iterable of three integers
    :param coo:  an indexible object that maps
    integers to x-y-z coordinates, e.g. a dictionary of tuples  
    """  
    x = []; y = []; 
    for p in range(3):
        vertex = triangle[ p ]
        x.append( coo[vertex][0] )
        y.append( coo[vertex][1] )
    return go.Scatter(x=x, y=y, fill="toself",)

def triangles__trace2d(triangles_as_rows=[], coo=[], single_trace=True):
    """
    Generates a Plotly `Mesh3d` trace.

    :param triangles_as_rows: an array or list of lists of
    integers, where each integer points to a row of 
    :param coo: an indexible object that maps
    integers to x-y-z coordinates, e.g. a dictionary of tuples  
    """

    # Plotly Mesh3d requires x, y, z to be formatted as lists, tuples, or arrays
    # for this reason, to avoid bundling a large point cloud into this trace, we
    # relabel the vertices incident to the simplices with labels {1, ..., n}, and
    # then generate coordinate array with only as many coordinates as we need
    if single_trace:
        x = []
        y = []
        for triangle in triangles_as_rows:    
            # generate a closed loop for each triangle    
            for p in range(4):
                vertex = triangle[ p % 3 ]
                x.append( coo[vertex][0] )
                y.append( coo[vertex][1] )
            # mark the break between triangles with None
            x.append(None)
            y.append(None)
        
        trace = go.Scatter(x=x, y=y, fill="toself") 
        # fig = go.Figure(go.Scatter(x=[0,1,2,0,None,3,3,5,5,3], y=[0,2,0,0,None,0.5,1.5,1.5,0.5,0.5], fill="toself"))    
        return trace 
    else:
        data = []
        for triangle in triangles_as_rows:    
            x = []
            y = []
            # generate a closed loop for each triangle    
            for p in range(4):
                vertex = triangle[ p % 3 ]
                x.append( coo[vertex][0] )
                y.append( coo[vertex][1] )
            # mark the break between triangles with None
            x.append(None)
            y.append(None)
            data.append( go.Scatter(x=x, y=y, fill="toself")  )
        return data         



def triangles__trace3d(triangles_as_rows=[], coo=[]):
    """
    Generates a Plotly `Mesh3d` trace.

    :param triangles_as_rows: is an array or list of lists of
    integers, where each integer points to a row of 
    :param coo: an indexible object that maps
    integers to x-y-z coordinates, e.g. a dictionary of tuples 
    :return trace: a single trace containing all triangles 

    **Note for developers** This method tends to be less flexible than repeated calls to `triangle__trace3d`,
    however the developers have kept it in the library thus (1) to preserve a well-tested function, and (2)
    to conserve at least one functional example of how to plot multiple disjoint edges in a single trace.    
    """

    print("Suggestion: consider using `triangle__trace3d` instead of `triangles__trace3d`, as the former often gives higher quality graphics.")
    # Plotly Mesh3d requires x, y, z to be formatted as lists, tuples, or arrays
    # for this reason, to avoid bundling a large point cloud into this trace, we
    # relabel the vertices incident to the simplices with labels {1, ..., n}, and
    # then generate coordinate array with only as many coordinates as we need
    triangles_as_rows = np.array(triangles_as_rows)
    nvl2ovl = np.unique( np.array(triangles_as_rows).flatten() )
    ovl2nvl = { nvl2ovl[nvl]: nvl for nvl in range(len(nvl2ovl))}
    triangles_as_rows_relabeled = [ [ ovl2nvl[x] for x in simplex] for simplex in triangles_as_rows ]
    triangles_as_rows_relabeled = np.array(triangles_as_rows_relabeled)

    coordinates_as_rows = np.array( [ coo[ v ] for v in nvl2ovl ] )

    trace   =   go.Mesh3d(
                    x=coordinates_as_rows[:,0],
                    y=coordinates_as_rows[:,1],
                    z=coordinates_as_rows[:,2],
                    i = [x[0] for x in triangles_as_rows_relabeled],
                    j = [x[1] for x in triangles_as_rows_relabeled],
                    k = [x[2] for x in triangles_as_rows_relabeled],
                )    
    return trace 



#   ==============================================================================
#   SHAPES 2D
#   ==============================================================================


def ball_2d( x, y, radius, npoints ):
    """
    Returns a trace for a ball of radius `radius` centered at `(x,y)`.  The perimiter
    of the ball is a closed piecewise linear curve with `npoints` vertices.
    """
    theta = np.linspace( 0, 2 * np.pi, npoints )    
    return go.Scatter(
        x   =   x + radius * np.cos(theta), 
        y   =   y + radius * np.sin(theta), 
        fill="toself"
    )




#   ==============================================================================
#   WIRE DIAGRAM 3D
#   ==============================================================================


def wire_sphere3d(x, y, z, radius, nlattitude, nlongitude):
    """
    Plot a wire mesh sphere

    Returns a list of traces which collectively make the wire mesh.

    - you can cause all this traces to toggle on/off together by calling `trace.update(legendgroup="your_legend_group_name")` on each trace;
    this is docuemented on the [plotly website](https://plotly.com/python/legend/) and you can
    find examples in the Exact gallery.
    - you can then remove all by one of the traces from the legend by calling `trace.update(showlegend=False)` on the traces you wish to hid
    """

    nlongitude = nlongitude + 1 # this corrects for the fact that we overlap 1 line

    # Define the phi and theta values for the sphere
    phi = np.linspace(0, np.pi, nlattitude)
    theta = np.linspace(0, 2*np.pi, nlongitude)

    # Compute the x, y, and z coordinates of the sphere points
    x = x + radius * np.outer(np.sin(phi), np.cos(theta)).ravel()
    y = y + radius * np.outer(np.sin(phi), np.sin(theta)).ravel()
    z = z + radius * np.outer(np.cos(phi), np.ones_like(theta)).ravel()

    # Define the indices of the points on each latitudinal line
    indices = [np.arange(i*nlongitude, (i+1)*nlongitude) for i in range(nlattitude)]

    # Define the indices of the points on each longitudinal line
    for i in range(nlongitude):
        indices.append(np.arange(i, nlattitude*nlongitude, nlongitude))

    #   Create the longitudinal lines
    data    =   [go.Scatter3d(x=x, y=y, z=z, mode='lines', line=dict(color='black', width=1))]

    # Add the latitudinal lines to the plot
    for index in indices:
        data.append(go.Scatter3d(x=x[index], y=y[index], z=z[index], mode='lines', line=dict(color='black', width=1)))    

    return data



#   ==============================================================================
#   SURFACE PLOTS
#   ==============================================================================



#   CUBE
#   --------------------------------------
def surface_cube(x0,y0,z0, width=1, anchor="center"):
    """
    A Plotly Surface trace for a 3-dimensional cube with side length `width`
    ```
    [x0,x0+width] x [y0,y0+width] x [z0,z0+width]
    ```
    :param x float: x coordinate of the anchor point
    :param y float: y coordinate of the anchor point
    :param z float: z coordinate of the anchor point        
    :param achor str: if `anchor="left" then the anchor point is the point on the cube with minimal x, y, and z coordinates.  If `anchor="center"` then the anchor is the centerpoint.
    :return trace, x, y, z: where `trace = go.Surface(z=z, x=x, y=y)`

    This is a convenience wrapper around the `surface_rectangle` function.
    """
    if anchor == "left":
        x1 = x0 + width
        y1 = y0 + width
        z1 = z0 + width
    elif anchor == "center":
        x0 = x0 - width/2
        y0 = y0 - width/2
        z0 = z0 - width/2
        x1 = x0 + width/2
        y1 = y0 + width/2
        z1 = z0 + width/2
    else:
        raise ValueError('The "anchor" keyword argument must be "left" or "center"')
        
    return surface_rectangle(x0,x1,y0,y1,z0,z1)


#   RECTANGLE
#   --------------------------------------
def surface_rectangle(x0,x1,y0,y1,z0,z1):
    """
    A Plotly Surface trace for a 3-dimensional rectangle of form [x0,x1] x [y0,y1] x [z0,z1].

    Returns one trace and three coordinate matrices: `go.Surface(z=z, x=x, y=y), x, y, z`.
    The coordinates tend to be useful for adjusting surface color.

    **Remark** This tends to produce cleaner results than a Plotly mesh plot.

    **How it works**

    The visual intuition is to image a cloth 3 units wide and 4 units long being
    folded over the surface of a cube.  The lefthand 3x3 units of cloth cover the top
    5 faces of the cube, and the righthand 1x3 units tuck under to cover the base.
    Some squares collapse down to lines. For reference, calling 

    ```
    _, x, y, z = rectangle_trace(0,1,2,3,4,5)
    print(x,"\n\n",y,"\n\n",z)
    ```

    will return

    ```
    [[0 0 1 1 0]
     [0 0 1 1 0]
     [0 0 1 1 0]
     [0 0 1 1 0]] 

    [[2 2 2 2 2]
     [2 2 2 2 2]
     [3 3 3 3 3]
     [3 3 3 3 3]] 

    [[4 4 4 4 4]
     [4 5 5 4 4]
     [4 5 5 4 4]
     [4 4 4 4 4]]
    ```
    """
    z = np.full((4,5), z0)
    z[1:3,1:3]= z1

    x = np.full((4,5), x0)
    x[:,2:4]=x1

    y = np.full((4,5), y0)
    y[2:,:]=y1

    return go.Surface(z=z, x=x, y=y), x, y, z,


#   SPHERE
#   --------------------------------------
def surface_sphere(x, y, z, radius, resolution=20):
    """Trace for a sphere"""
    u, v = np.mgrid[0:2*np.pi:resolution*2j, 0:np.pi:resolution*1j]
    x = radius * np.cos(u)*np.sin(v) + x
    y = radius * np.sin(u)*np.sin(v) + y
    z = radius * np.cos(v) + z
    trace = go.Surface(x=x,y=y,z=z), x, y, z
    return trace

    


#   ==============================================================================
#   MESH PLOTS
#   ==============================================================================



#   RECTANGLE TRACE (MESH - WORSE RESULTS)
#   --------------------------------------
def mesh_rectangle(x0,x1,y0,y1,z0,z1):
    """
    A Plotly Mesh3d trace for a rectangle of form [x0,x1] x [y0,y1] x [z0,z1].
    """
    print("This method tends to generate poor shadowing and seams at non-right angles.  Consider `rectangle_trace_3d` as an alternative.")

    x = np.array([0, 0, 1, 1, 0, 0, 1, 1]) * (x1-x0) + x0
    y = np.array([0, 1, 1, 0, 0, 1, 1, 0]) * (y1-y0) + y0
    z = np.array([0, 0, 0, 0, 1, 1, 1, 1]) * (z1-z0) + z0   

    return go.Mesh3d(
        # 8 vertices of a cube
        x=x,
        y=y,
        z=z,
        colorbar_title='z',
        colorscale=[[0, 'gold'],
                    [0.5, 'mediumturquoise'],
                    [1, 'magenta']],
        # Intensity of each vertex, which will be interpolated and color-coded
        intensity = np.linspace(0, 1, 8, endpoint=True),
        # i, j and k give the vertices of triangles
        i = [7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2],
        j = [3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3],
        k = [0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6],
        name='y',
        showscale=True
    ) 



def surface_octahedron():
    """
    Returns a trace composed of 8 triangles forming an octahedron
    """
    coo = np.array( # coordinate oracle
        [
            # first four columns "walk around the equator along adjacent vertices"
            # the final two columns represent the north/south poles
            [ -1,  0, 1,  0,  0,  0 ], # x
            [  0, -1, 0,  1,  0,  0 ], # y
            [  0,  0, 0,  0, -1,  1 ], # z
        ]
    ).T

    edges = [ [0,1], [1,2], [2,3], [3,0] ]
    triangles = [ edge + [pole] for edge in edges for pole in [4,5] ]

    data = []
    for triangle in triangles:
        trace = oat_python.plot.triangle__trace3d( triangle, coo=coo)
        trace.update(name=f"{triangle}", showlegend=True, color="red", opacity=0.5)
        data.append(trace)
    return data






#   ==============================================================================
#   UNDER REVIEW FOR DEPRECATION
#   ==============================================================================


#   !!! BEFORE DEPRECATING, CONSIDER THE RELATIONSHIP BETWEEN THIS FUNCTION AND MDS; IT MAY BE THE ONLY PLACE WE USE OUR MDS FUNCTIONS
def plot_simplices( 
        simplices=[],     
        coo=None, 
        simplex_labels=None, 
        vertex_labels=None,
        ovl2intensity=None,
        showcloud=False
    ):
    """
    Plot a simplicial complex in interactive 3d, with Plotly.

    The user specifies a list of simplices (represented as a list of lists of integers) and x-y-z coordinates for each vertex.

    :param coo either None (generate coordinates automatically using MDS) or numpy.array coo: array of size p x 3, where p is the number of vertices; row k represents the x,y,z coordinates of the kth vertex
    :param iterable of lists: each list represents a simplex
    :param simplex_labels (optional): list of labels for each simplex; a label can be anything that the `str` function can convert to a string; if this argument is not specified then simplices will be labeled as lists of vertices
    :param vertex_labels (optional): list of labels for each vertex; a label can be anything that the `str` function can convert to a string; if this argument is not specified then vertices wil be number 0 through p-1, where k corresponds to the kth row of the matrix coo
    :param ovl2intensity (optional): map sending old vertex labels to intensity
    """

    print("`plot_simplices` is under review for deprecation, as it often fails to give sufficiently fine-grained control.")
    simplices = np.array(simplices)
    nvl2ovl = np.unique( np.array(simplices).flatten() )
    ovl2nvl = { nvl2ovl[nvl]: nvl for nvl in range(len(nvl2ovl))}
    simplices_relabled = [ [ ovl2nvl[x] for x in simplex] for simplex in simplices ]
    simplices_relabled = np.array(simplices_relabled)

    if ovl2intensity is None:
        ovl2intensity = np.zeros(np.max(nvl2ovl)+1)

    # supply coordinates, if not provided
    if coo is None:
        coo = hop_mds_from_simplices(simplices)

    coordinates_as_rows = np.array( [coo[k] for k in nvl2ovl] )        
    

    #   a function to generate simplex labels for the legend
    def make_simplex_name( simplex_num ):
        if simplex_labels is None:
            simplex = simplices_relabled[simplex_num]
            if vertex_labels is None:
                return str( simplex )
            else:
                return ','.join([ str(vertex_labels[p]) for p in simplex])
        else:
            return str(simplex_labels[simplex_num])

    # create a list of plots
    data    =   []

    # add a scatter for the point cloud
    if vertex_labels is None:
        vertex_labels   =   nvl2ovl
    if showcloud:
        trace = go.Scatter3d(
            x=coordinates_as_rows[:,0], 
            y=coordinates_as_rows[:,1],
            z=coordinates_as_rows[:,2],
            mode='markers+text',
            text=vertex_labels,
            name="Vertices"
            )
        data.append(trace)

    # helper function to format the list of polytopes that we'll pass to the meshgrid function to generate plots of 2-simplices
    def cardinality3subsets_as_rows(simplex):
        """
        returns np array whos rows are all cardinaltiy-three subsets of the given set
        """
        return  np.array( [ combo for combo in itertools.combinations( simplex, 3 ) ] )


    # generate plots for all simplices
    for simplex_num, simplex in enumerate(simplices_relabled):
        if len(simplex) == 0:
            vertex  =   simplex[0]
            trace   =   go.Scatter3d(
                            x=[], 
                            y=[],
                            z=[],
                            mode='markers',
                            name= make_simplex_name(simplex_num),
                        )
        elif len(simplex) == 1:
            vertex  =   simplex[0]
            trace   =   go.Scatter3d(
                            x=[coordinates_as_rows[vertex,0]], 
                            y=[coordinates_as_rows[vertex,1]],
                            z=[coordinates_as_rows[vertex,2]],
                            mode='markers',
                            name= make_simplex_name(simplex_num),
                        )
        elif len(simplex) == 2:
            vertex  =   simplex[0]
            trace   =   go.Scatter3d(
                            x=coordinates_as_rows[simplex,0], 
                            y=coordinates_as_rows[simplex,1],
                            z=coordinates_as_rows[simplex,2],
                            mode='lines',
                            name= make_simplex_name(simplex_num),
                        )                        
        elif len(simplex) >= 3:
            coords  =   coordinates_as_rows[simplex,:]            
            card3subsets    =   cardinality3subsets_as_rows( range(len(simplex)) ) # we do this because we've re-indexed the array of coordinates
            trace   =   go.Mesh3d(
                            x=coords[:,[0]].flatten(),
                            y=coords[:,[1]].flatten(),
                            z=coords[:,[2]].flatten(),
                            # colorbar_title='z',
                            # colorscale=[[0, 'gold'],
                            #             [0.5, 'mediumturquoise'],
                            #             [1, 'magenta']],
                            # # Intensity of each vertex, which will be interpolated and color-coded
                            # intensity=[0, 0.33, 0.66, 1],
                            # i, j and k give the vertices of triangles
                            # here we represent the 4 triangles of the tetrahedron surface
                            i= card3subsets[:,[0]].flatten(),
                            j= card3subsets[:,[1]].flatten(),
                            k= card3subsets[:,[2]].flatten(),
                            opacity=0.5,
                            intensity= [ ovl2intensity[ nvl2ovl[ x ] ] for x in simplex ],
                            # cmin=0.0,
                            # cmax=1.0,
                            name= make_simplex_name(simplex_num),
                            showlegend=True,

                        )                     
        data.append(trace)
    
    return data