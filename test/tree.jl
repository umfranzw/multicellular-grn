using TreeVisMod
using Gtk

file_path = "/home/wayne/test.png"

code = """
digraph G {
  1 [label="1"]
  1 -> 2
  1 -> 3
  2 [label="2"]
  3 [label="3"]
}
"""

data = TreeVisMod.gen_graph(code)
f = open(file_path, "w")
write(f, data)
close(f)

image = GtkImage(file_path)
w = GtkWindow(image, "title")
show(image)
