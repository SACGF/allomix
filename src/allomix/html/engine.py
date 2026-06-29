"""Jinja2 environment for the HTML report.

Builds the template environment and registers the formatting helpers from
``allomix.html.format`` as Jinja filters/globals, so templates do no analysis or
number formatting of their own: they only place already-formatted, HTML-safe
values. Autoescaping is on, so every value interpolated with ``{{ ... }}`` is
escaped unless it is explicitly marked safe (only the status badge, which emits a
fixed span, is).

The loader is a ``ChoiceLoader``: an optional user template directory is searched
first, then the built-in templates that ship with the package. A lab can drop a
single file (for example ``styles.css`` to restyle, or ``report.html`` to
restructure) into its own directory and override just that file, falling back to
the built-ins for everything else (the ``--template`` CLI flag).
"""

from pathlib import Path

from jinja2 import ChoiceLoader, Environment, FileSystemLoader, PackageLoader, select_autoescape
from markupsafe import Markup

from allomix.html import format as fmt

# Templates (and the CSS/JS assets) ship as package data under this directory.
_PACKAGE = "allomix"
_TEMPLATE_PACKAGE_PATH = "html/templates"


def make_environment(template_dir: str | Path | None = None) -> Environment:
    """Build the report's Jinja2 environment.

    Args:
        template_dir: Optional directory of user template overrides, searched
            ahead of the built-in templates. Any file present there (a full
            ``report.html``, or just ``styles.css``) overrides the packaged one;
            anything absent falls back to the built-in.

    Returns:
        A configured ``Environment`` with the formatting filters/globals
        registered.
    """
    loaders = []
    if template_dir is not None:
        loaders.append(FileSystemLoader(str(template_dir)))
    loaders.append(PackageLoader(_PACKAGE, _TEMPLATE_PACKAGE_PATH))

    env = Environment(
        loader=ChoiceLoader(loaders),
        autoescape=select_autoescape(["html", "xml"]),
        trim_blocks=True,
        lstrip_blocks=True,
    )

    # Number / percentage / p-value formatters, as filters: ``{{ x | pct(3) }}``.
    env.filters["num"] = fmt.num
    env.filters["count"] = fmt.count
    env.filters["pct"] = fmt.pct
    env.filters["pct_points"] = fmt.pct_points
    env.filters["pval"] = fmt.pval
    env.filters["pct_value"] = _pct_value

    # Helpers taking more than the piped value, as globals: ``{{ ci(lo, hi) }}``.
    env.globals["ci"] = fmt.ci
    env.globals["pct_value"] = _pct_value
    env.globals["NA"] = fmt.NA
    env.globals["badge"] = _badge

    return env


def load_asset(env: Environment, name: str) -> Markup:
    """Return a static asset's raw text (CSS/JS), honouring template overrides.

    The asset is read through the environment loader (so a user ``template_dir``
    can override it) and returned as ``Markup`` so it is inlined verbatim rather
    than HTML-escaped. It is not parsed as a template, so user CSS/JS may contain
    any characters, including Jinja delimiters.

    Args:
        env: The environment from ``make_environment``.
        name: Asset template name, for example ``"styles.css"``.

    Returns:
        The asset text, marked safe for inlining.
    """
    source, _, _ = env.loader.get_source(env, name)
    return Markup(source)


def _badge(status: str, large: bool = False) -> Markup:
    """Status badge as safe HTML (the only caller-independent HTML helper)."""
    return Markup(fmt.badge(status, large=large))


def _pct_value(v: float | None, dp: int = 2) -> str:
    """Format an already-percent value (0-100), or NA when not finite."""
    if v is None or not (v == v and v not in (float("inf"), float("-inf"))):
        return fmt.NA
    return f"{v:.{dp}f}%"
