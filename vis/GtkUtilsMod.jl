module GtkUtilsMod

using Gtk
using Gtk.GLib
using Plots

#note: Gtk v0.18.0 (ibis) Gtk.libgdk_pixbuf
#Gtk v1.1.2 (laptop) Gtk.libgdkpixbuf

function plot_to_pixbuf(plot::Plots.Plot)
    buf = IOBuffer()
    show(buf, MIME("image/png"), plot)

    pixbuf_from_data(buf.data)
end

function pixbuf_from_data(data::Array{UInt8, 1})
    pixbuf_loader = ccall((:gdk_pixbuf_loader_new_with_type, Gtk.libgdkpixbuf),
                          Ptr{GObject},
                          (Ptr{UInt8}, Ptr{Ptr{GError}}),
                          GLib.bytestring("png"), C_NULL
                          )

    result = ccall((:gdk_pixbuf_loader_write, Gtk.libgdkpixbuf),
                   Bool,
                   (Ptr{GObject}, Ptr{UInt8}, UInt64, Ptr{Ptr{GError}}),
                   pixbuf_loader, data, UInt64(length(data)), C_NULL
                   )

    pixbuf = ccall((:gdk_pixbuf_loader_get_pixbuf, Gtk.libgdkpixbuf),
                   Ptr{GObject},
                   (Ptr{GObject},),
                   pixbuf_loader
                   )

    GdkPixbufLeaf(pixbuf)
end

end
