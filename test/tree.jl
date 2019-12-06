using Gtk

ts = GtkTreeStore(String)
iter1 = push!(ts,("one",))
iter2 = push!(ts,("two",),iter1)
iter3 = push!(ts,("three",),iter2)
tv = GtkTreeView(GtkTreeModel(ts))
r1 = GtkCellRendererText()
c1 = GtkTreeViewColumn("A", r1, Dict([("text",0)]))
push!(tv,c1)
win = GtkWindow(tv, "Tree View")
showall(win)

iter = Gtk.iter_from_index(ts, [1])
ts[iter,1] = "ONE"
