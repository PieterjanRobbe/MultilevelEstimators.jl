#******************************************************************************
# WARNING: do not change the lines below if you want support from the staff
#******************************************************************************

if [ -f /usr/skel/new/profile.called ]; then
	CS_PROFILE=yes
	. /usr/skel/new/profile.called
else
	echo "+----------------------------------------------------------------------+"
	echo "| WARNING : No system dependent defaults found, using minimal defaults |"
	echo "+----------------------------------------------------------------------+"
	PATH=/usr/local/bin:/usr/bin:/bin
	export PATH
	MANPATH=/usr/local/man:/usr/man
	export MANPATH
fi
#******************************************************************************
# WARNING: do not change the lines above if you want support from the staff
#******************************************************************************

#--------------------------------------------------------------------
# This file was last modified on mar 20th, 1997 by a member of staff.
# You can check for a new version in the /usr/skel/new directory.
#--------------------------------------------------------------------

#****************************#
# do your private stuff here #
#****************************#
#
# The "systeemgroep" strongly advises against 
# indiscriminately copying of major parts of this file
# between users, because the method is error prone, and
# the resulting problems are often difficult to trace.
# If you do this, you better know what you are doing!
#
# However, if you think that your stuff might be useful for many other
# users, activate your representative in the systems group for inclusion
# of your enhancements/features in the global files.
# 
# Also, please do report all problems you might experience with this
# setup.  We cannot give good support for this kind of files without
# your cooperation.
#
# If you mess up this file too much, you can get a fresh copy
# with following command : cp /usr/skel/new/bash_profile ${HOME}/.bash_profile

MATLAB_HOME=/cw/matlab-2018b-64
export MATLAB_HOME

# From here on you have a number of useful variables set
# They include :
#  PATH
#  MANPATH
#  USER
#  ...
# The user changeable options are stored in
# your home-directory in the file ".options.sh"
# These options control the contents of above mentioned variables.

# If your .options.sh is lost or is outdated, you can copy
# a fresh one with following command:
# cp /usr/skel/new/options.sh $HOME/.options.sh

