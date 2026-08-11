import unittest
from pathlib import Path
from unittest.mock import patch

from gui import theme


ROOT = Path(__file__).resolve().parents[1]


class RememberedTabTests(unittest.TestCase):
    def test_widget_callback_records_owning_tab(self):
        state = {}
        with patch.object(theme.st, "session_state", state):
            theme.remember_tab("active_tab", "Plot Styling")
        self.assertEqual(state["active_tab"], "Plot Styling")

    def test_tabs_reopen_at_remembered_value(self):
        state = {"active_tab": "Second"}
        with patch.object(theme.st, "session_state", state):
            with patch.object(
                theme.st, "tabs", return_value=("one", "two")
            ) as tabs:
                rendered = theme.remembered_tabs(
                    ["First", "Second"], "active_tab"
                )
        self.assertEqual(rendered, ("one", "two"))
        tabs.assert_called_once_with(["First", "Second"], default="Second")

    def test_stale_tab_value_falls_back_safely(self):
        state = {"active_tab": "Removed tab"}
        with patch.object(theme.st, "session_state", state):
            with patch.object(theme.st, "tabs", return_value=()) as tabs:
                theme.remembered_tabs(["First", "Second"], "active_tab")
        self.assertEqual(state["active_tab"], "First")
        tabs.assert_called_once_with(["First", "Second"], default="First")

    def test_pages_use_shared_remembered_tab_helper(self):
        direct_users = []
        for path in sorted((ROOT / "gui").glob("page_*.py")):
            if "st.tabs(" in path.read_text(encoding="utf-8"):
                direct_users.append(path.name)
        self.assertEqual(direct_users, [])


if __name__ == "__main__":
    unittest.main()
