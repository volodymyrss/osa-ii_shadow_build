#ifndef __SHADOW_BUILD_H__
#define __SHADOW_BUILD_H__

#include "ii_shadow_build_f.h"

{% if mode == "delivery" %}
#define COMPONENT_NAME           "{{component_name}}"
#define COMPONENT_VERSION        "{{component_version}}"
{% else %}
#ifndef GIT_VERSION
#define GIT_VERSION ""
#endif
#define COMPONENT_NAME           "{{component_name}}"
#define COMPONENT_VERSION        "{{component_version}}-"GIT_VERSION
{% endif %}


#endif
