import RunMod
using TrackerMod

run = RunMod.get_first_run()
TrackerMod.create_tracker(run)
tracker = TrackerMod.get_tracker()
