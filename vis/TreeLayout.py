class TreeLayout():
    @staticmethod
    def buchheim(draw_tree):
        dt = TreeLayout.firstwalk(draw_tree)
        min = TreeLayout.second_walk(dt)
        if min < 0:
            TreeLayout.third_walk(dt, -min)
        return dt

    @staticmethod
    def third_walk(tree, n):
        tree.x += n
        for c in tree.children:
            TreeLayout.third_walk(c, n)

    @staticmethod
    def firstwalk(v, distance=1.):
        if len(v.children) == 0:
            if v.lmost_sibling:
                v.x = v.lbrother().x + distance
            else:
                v.x = 0.
        else:
            default_ancestor = v.children[0]
            for w in v.children:
                TreeLayout.firstwalk(w)
                default_ancestor = TreeLayout.apportion(w, default_ancestor, distance)
            #print("finished v =", v.tree, "children")
            TreeLayout.execute_shifts(v)

            midpoint = (v.children[0].x + v.children[-1].x) / 2

            ell = v.children[0]
            arr = v.children[-1]
            w = v.lbrother()
            if w:
                v.x = w.x + distance
                v.mod = v.x - midpoint
            else:
                v.x = midpoint
        return v

    @staticmethod
    def apportion(v, default_ancestor, distance):
        w = v.lbrother()
        if w is not None:
            #in buchheim notation:
            #i == inner; o == outer; r == right; l == left; r = +; l = -
            vir = vor = v
            vil = w
            vol = v.lmost_sibling
            sir = sor = v.mod
            sil = vil.mod
            sol = vol.mod
            while vil.right() and vir.left():
                vil = vil.right()
                vir = vir.left()
                vol = vol.left()
                vor = vor.right()
                vor.ancestor = v
                shift = (vil.x + sil) - (vir.x + sir) + distance
                if shift > 0:
                    TreeLayout.move_subtree(TreeLayout.ancestor(vil, v, default_ancestor), v, shift)
                    sir = sir + shift
                    sor = sor + shift
                sil += vil.mod
                sir += vir.mod
                sol += vol.mod
                sor += vor.mod
            if vil.right() and not vor.right():
                vor.thread = vil.right()
                vor.mod += sil - sor
            else:
                if vir.left() and not vol.left():
                    vol.thread = vir.left()
                    vol.mod += sir - sol
                default_ancestor = v
        return default_ancestor

    @staticmethod
    def move_subtree(wl, wr, shift):
        subtrees = wr.number - wl.number
        #print(wl.tree, "is conflicted with", wr.tree, 'moving', subtrees, 'shift', shift)
        #print wl, wr, wr.number, wl.number, shift, subtrees, shift/subtrees
        wr.change -= shift / subtrees
        wr.shift += shift
        wl.change += shift / subtrees
        wr.x += shift
        wr.mod += shift

    @staticmethod
    def execute_shifts(v):
        shift = change = 0
        for w in v.children[::-1]:
            #print("shift:", w, shift, w.change)
            w.x += shift
            w.mod += shift
            change += w.change
            shift += w.shift + change

    @staticmethod
    def ancestor(vil, v, default_ancestor):
        #the relevant text is at the bottom of page 7 of
        #"Improving Walker's Algorithm to Run in Linear Time" by Buchheim et al, (2002)
        #http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.16.8757&rep=rep1&type=pdf
        if vil.ancestor in v.parent.children:
            return vil.ancestor
        else:
            return default_ancestor

    @staticmethod
    def second_walk(v, m=0, depth=0, min=None):
        v.x += m
        v.y = depth

        if min is None or v.x < min:
            min = v.x

        for w in v.children:
            min = TreeLayout.second_walk(w, m + v.mod, depth+1, min)

        return min
