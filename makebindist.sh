#!/bin/tcsh

setenv NAME miRanda-`cat VERSION`-`./config.guess | sed 's/^\([^-]*\)-\([^-]*\)-\(.*\)$/\1/'`-`./config.guess | sed 's/^\([^-]*\)-\([^-]*\)-\(.*\)$/\3/'`
echo ${NAME}

./configure --prefix ${PWD}/${NAME}
make install
cp LICENSE README COPYING AUTHORS NEWS THANKS ${NAME}
cp -R examples ${NAME}
cp man/miranda.* ${NAME}/man
tar cvf miRanda-`cat VERSION`-`./config.guess | sed 's/^\([^-]*\)-\([^-]*\)-\(.*\)$/\1/'`-`./config.guess | sed 's/^\([^-]*\)-\([^-]*\)-\(.*\)$/\3/'`.tar ${NAME}
gzip miRanda-`cat VERSION`-`./config.guess | sed 's/^\([^-]*\)-\([^-]*\)-\(.*\)$/\1/'`-`./config.guess | sed 's/^\([^-]*\)-\([^-]*\)-\(.*\)$/\3/'`.tar
rm -rf ${NAME}

