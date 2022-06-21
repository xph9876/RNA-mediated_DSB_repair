import plot_graph
import plotly.graph_objects as go
import plotly.subplots as ps

# tick_pos = [
#   {'pos': 1, 'text': '1'},
#   {'pos': 2, 'text': '2'},
#   {'pos': 3, 'text': '3'},
#   {'pos': 4, 'text': '4'},
#   {'pos': 5, 'text': '5'},
# ]
# fig = go.Figure()
# fig.add_shape(
#   x0 = 0,
#   x1 = 0,
#   y0 = 0,
#   y1 = 5,
# )
# for tick in tick_pos:
#   fig.add_shape(
#     type = 'line',
#     x0 = 0,
#     x1 = 0.75,
#     y0 = tick['pos'],
#     y1 = tick['pos'],
#   )
#   fig.add_annotation(
#     x = 1,
#     y = tick['pos'],
#     text = tick['text'],
#     showarrow = False,
#   )
# fig.show()

fig = ps.make_subplots(1, 1)
plot_graph.make_universal_layout_legend(fig, 0, 1, 1)
fig.update_yaxes(range = [0, 20])
fig.show()