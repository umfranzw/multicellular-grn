using Gtk
using Gadfly
using Compose

p = plot(x=rand(40, 1), y=rand(40, 1))
co = render(p)
c = GtkCanvas(400, 300)
win = GtkWindow(c, "data win")
show(c)

Gtk.draw(c) do widget
    Compose.draw(CAIROSURFACE(c.back), co)
end

condition = Condition()
endit(w) = notify(condition)
signal_connect(endit, win, :destroy)
showall(win)
wait(condition)
