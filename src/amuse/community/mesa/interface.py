import amuse
if amuse.config.community.mesa_version == "2208":
    from amuse.community.mesa_2208.interface import *
elif amuse.config.community.mesa_version == "new":
    from amuse.community.mesa_new.interface import *
