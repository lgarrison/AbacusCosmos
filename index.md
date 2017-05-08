---
layout: default
---

{% for item in site.data.navigation.pages %}
<li class="{% if item.url == page.url %}active{% endif %}">
  <a href="{{ item.url }}">{{ item.title }}</a>
</li>
{% endfor %}
