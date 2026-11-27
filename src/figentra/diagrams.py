"""Conceptual scientific diagrams with editable vector output."""
from math import cos, sin, pi
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch
from .style import style_context

__all__ = ["venn_diagram", "research_workflow", "problem_tree",
           "conceptual_framework", "network_diagram", "biological_pathway",
           "mechanism_diagram"]

def _canvas(size_mm, background):
    with style_context():
        fig, ax = plt.subplots(figsize=tuple(v / 25.4 for v in size_mm))
    fig.patch.set_facecolor(background)
    ax.set_facecolor(background)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    fig.subplots_adjust(left=.025, right=.975, bottom=.025, top=.975)
    return fig, ax

def venn_diagram(labels, *, center="", colors=("#3979B9", "#C65F4A", "#38875C"),
                 size_mm=(160, 140), font_size=12, font_weight="bold",
                 line_width=2, alpha=.18, background="white", text_color="black"):
    """Draw a conceptual three-circle overlap, not an area-proportional Venn."""
    if len(labels) != 3 or len(colors) != 3:
        raise ValueError("Exactly three labels and colors are required")
    fig, ax = _canvas(size_mm, background)
    ax.set_aspect("equal")
    for (x,y), (tx,ty), label, color in zip(
        [( .33,.69),(.67,.69),(.50,.39)],
        [(.24,.76),(.76,.76),(.50,.25)], labels, colors):
        ax.add_patch(Circle((x,y), .30, facecolor=color, alpha=alpha, edgecolor="none"))
        ax.add_patch(Circle((x,y), .30, facecolor="none", edgecolor=color, linewidth=line_width))
        ax.text(tx,ty,label,ha="center",va="center",fontsize=font_size,
                fontweight=font_weight,color=text_color)
    if center:
        ax.text(.5,.59,center,ha="center",va="center",fontsize=font_size,
                fontweight=font_weight,color="white",
                bbox=dict(boxstyle="round,pad=.4",facecolor="#A63838",edgecolor="none"))
    return fig, ax

def _node(ax, x, y, label, w, h, color, font_size, font_weight, line_width, text_color):
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle="round,pad=0.008",
                              facecolor=color,edgecolor="#526477",linewidth=line_width,zorder=2))
    ax.text(x,y,label,ha="center",va="center",fontsize=font_size,
            fontweight=font_weight,color=text_color,zorder=3)

def _arrow(ax, a, b, line_width):
    ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="->",
                color="#526477",lw=line_width),zorder=1)

def research_workflow(steps, *, size_mm=(150,180), font_size=12,
                      font_weight="bold", line_width=2, color="#E5EFF8",
                      background="white", text_color="black"):
    """Draw a vertical sequence. Use line breaks for long labels."""
    steps = list(steps)
    if not steps:
        raise ValueError("At least one step is required")
    fig, ax = _canvas(size_mm, background)
    n = len(steps)
    h = min(.13,.65/n)
    ys = [1-(i+.5)/n for i in range(n)]
    for i,(y,label) in enumerate(zip(ys,steps)):
        _node(ax,.5,y,label,.8,h,color,font_size,font_weight,line_width,text_color)
        if i:
            _arrow(ax,(.5,ys[i-1]-h/2),(.5,y+h/2),line_width)
    return fig, ax

def problem_tree(problem, causes, consequences, *, size_mm=(200,140),
                 font_size=11, font_weight="bold", line_width=2,
                 background="white", text_color="black",
                 colors=("#E5F1E8","#E5EFF8","#FCE8E5")):
    """Draw causes -> core problem -> consequences."""
    causes, consequences = list(causes), list(consequences)
    if not causes or not consequences:
        raise ValueError("Provide causes and consequences")
    fig, ax = _canvas(size_mm, background)
    _node(ax,.5,.5,problem,.60,.14,colors[1],font_size,font_weight,line_width,text_color)
    for labels,y,color in [(causes,.15,colors[0]),(consequences,.85,colors[2])]:
        for i,label in enumerate(labels):
            x=(i+.5)/len(labels)
            _node(ax,x,y,label,.88/len(labels),.15,color,font_size,font_weight,line_width,text_color)
            a,b=((x,y+.075),(.5,.43)) if y<.5 else ((.5,.57),(x,y-.075))
            _arrow(ax,a,b,line_width)
    return fig, ax

def network_diagram(nodes, edges, *, size_mm=(180,150), font_size=11,
                    font_weight="bold", line_width=2, color="#E5EFF8",
                    background="white", text_color="black", directed=False):
    """Draw a circular network from unique node labels and edge pairs."""
    nodes=list(nodes)
    if not nodes or len(set(nodes))!=len(nodes):
        raise ValueError("Nodes must be nonempty and unique")
    pos={n:(.5+.36*cos(2*pi*i/len(nodes)),.5+.36*sin(2*pi*i/len(nodes))) for i,n in enumerate(nodes)}
    edges=list(edges)
    if any(a not in pos or b not in pos or a==b for a,b in edges):
        raise ValueError("Edges must connect distinct known nodes")
    fig, ax = _canvas(size_mm, background)
    for a,b in edges:
        if directed:
            ax.annotate("",xy=pos[b],xytext=pos[a],arrowprops=dict(arrowstyle="->",color="#526477",lw=line_width,shrinkA=24,shrinkB=24))
        else:
            ax.plot([pos[a][0],pos[b][0]],[pos[a][1],pos[b][1]],color="#526477",lw=line_width)
    for n,(x,y) in pos.items():
        _node(ax,x,y,n,.22,.09,color,font_size,font_weight,line_width,text_color)
    return fig, ax

def conceptual_framework(nodes, edges, **kwargs):
    """Directed conceptual relationships in a circular layout."""
    return network_diagram(nodes, edges, directed=True, **kwargs)

def biological_pathway(steps, **kwargs):
    """Draw a linear biological pathway."""
    return research_workflow(steps, **kwargs)

def mechanism_diagram(steps, **kwargs):
    """Draw a linear conceptual mechanism (not molecular structures)."""
    return research_workflow(steps, **kwargs)
