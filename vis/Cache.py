class Cache():
    def __init__(self, size):
        self.size = size
        self.insert_order = []
        self.data = {}

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, value):
        if key not in self.data:
            if len(self.data) + 1 > self.size:
                #remove the oldest item (key at the back of insert_order)
                oldest_key = self.insert_order.pop()
                del self.data[oldest_key]
            self.insert_order.insert(0, key)
        self.data[key] = value

    def keys(self):
        return self.data.keys()

    def __contains__(self, key):
        return key in self.data
