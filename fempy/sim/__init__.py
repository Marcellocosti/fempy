<<<<<<< HEAD
import os
# If you want to import alice3 sw add `export FEMPY_ALICE3=true` to your .bashrc
=======
# If you want to import alice3 sw add `export FEMPY_ALICE3=true` to your .bashrc
import os
>>>>>>> bef2228 (add alice3 software options)
if os.getenv("FEMPY_ALICE3", 'false').lower() in ('true', '1'):
    from . import alice3
