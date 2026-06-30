"""Optional report metadata for the HTML chimerism report.

These dataclasses carry the recipient / donor / transplant context shown in the
report header band. Everything is optional: the report must render cleanly when
metadata is absent and must never block on missing demographics. No values here
feed the analysis; they are presentational only.

They live in their own module (not ``report.py``) so both ``report.py`` and the
``allomix.report.html`` renderer can import them without an import cycle. ``report.py``
re-exports ``ReportMeta`` and ``DonorMeta`` as the public surface.
"""

from dataclasses import dataclass, field


@dataclass
class DonorMeta:
    """Identification of one donor for the report header.

    Attributes:
        donor_id: Donor identifier, as the lab labels it.
        relationship: Declared relationship to the recipient (for example
            "unrelated", "first-degree", "sibling"). Free text; shown verbatim.
            This is the human-facing label, separate from the
            ``--expected-relatedness`` value that drives the QC check.
    """

    donor_id: str | None = None
    relationship: str | None = None

    def to_dict(self) -> dict:
        """Serialise to a plain JSON-safe dict."""
        return {"donor_id": self.donor_id, "relationship": self.relationship}

    @classmethod
    def from_dict(cls, d: dict) -> "DonorMeta":
        """Rebuild from a ``to_dict`` mapping (unknown keys ignored)."""
        return cls(donor_id=d.get("donor_id"), relationship=d.get("relationship"))


@dataclass
class ReportMeta:
    """Optional recipient / transplant context for the report header band.

    All fields are optional. Dates are passed as already-formatted strings (ISO
    or otherwise); the report does not parse them. Days-post-transplant is
    derived only when both a transplant date and a sample date are supplied, and
    the derivation happens in the render layer (presentational only).

    Attributes:
        recipient_id: Recipient identifier.
        recipient_name: Optional recipient display name.
        sex: Optional recipient sex, shown verbatim.
        dob: Optional date of birth, as a preformatted string.
        transplant_type: Transplant type label (default "HSCT"); allows others.
        transplant_date: Optional transplant date, as a preformatted string.
        donors: Per-donor identification and declared relationship.
        sample_dates: Mapping of sample name to a preformatted collection-date
            string, used for the header rows and the timeline x-axis.
    """

    recipient_id: str | None = None
    recipient_name: str | None = None
    sex: str | None = None
    dob: str | None = None
    transplant_type: str = "HSCT"
    transplant_date: str | None = None
    donors: list[DonorMeta] = field(default_factory=list)
    sample_dates: dict[str, str] = field(default_factory=dict)

    def to_dict(self) -> dict:
        """Serialise to a plain JSON-safe dict for the report envelope."""
        return {
            "recipient_id": self.recipient_id,
            "recipient_name": self.recipient_name,
            "sex": self.sex,
            "dob": self.dob,
            "transplant_type": self.transplant_type,
            "transplant_date": self.transplant_date,
            "donors": [d.to_dict() for d in self.donors],
            "sample_dates": dict(self.sample_dates),
        }

    @classmethod
    def from_dict(cls, d: dict | None) -> "ReportMeta":
        """Rebuild from a ``to_dict`` mapping; an empty/None mapping yields the
        all-default metadata so a report renders cleanly without it."""
        if not d:
            return cls()
        return cls(
            recipient_id=d.get("recipient_id"),
            recipient_name=d.get("recipient_name"),
            sex=d.get("sex"),
            dob=d.get("dob"),
            transplant_type=d.get("transplant_type", "HSCT"),
            transplant_date=d.get("transplant_date"),
            donors=[DonorMeta.from_dict(x) for x in d.get("donors", [])],
            sample_dates=dict(d.get("sample_dates", {})),
        )
