<<<<<<< HEAD
<<<<<<< HEAD
import os
# If you want to import alice3 sw add `export FEMPY_ALICE3=true` to your .bashrc
import os
if os.getenv("FEMPY_ALICE3", 'false').lower() in ('true', '1'):
    from . import alice3
=======
#from . import alice3
>>>>>>> 4b5d47b (comment alice3 config)
=======
# If you want to import alice3 sw add `export FEMPY_ALICE3=true` to your .bashrc
if os.getenv("FEMPY_ALICE3", 'false').lower() in ('true', '1'):
    from . import alice3
>>>>>>> bef2228 (add alice3 software options)
