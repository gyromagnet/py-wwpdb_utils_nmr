"""A rejected source must still have a schema-valid diagnostic representation."""
import json
import unittest
from pathlib import Path

from jsonschema import Draft7Validator, ValidationError


class ReportSchemaRegressionTests(unittest.TestCase):

    def test_unidentified_content_is_explicitly_unknown(self):
        path = Path(__file__).parent.parent / 'tests-nmr' / 'json-schema' / 'nmr-data-processing-report-schema-v5.json'
        schema = json.loads(path.read_text(encoding='utf-8'))
        Draft7Validator.check_schema(schema)
        validator = Draft7Validator(schema['definitions']['input_source']['properties']['content_subtype'])
        validator.validate(None)
        validator.validate({'dist_restraint': 1})
        # Unknown is distinct from a known set of zero observations.
        validator.validate({})
        for invalid in ([], 'unknown', 0, {'dist_restraint': -1}):
            with self.subTest(invalid=invalid), self.assertRaises(ValidationError):
                validator.validate(invalid)

    def test_output_statistics_unknown_metadata_and_unavailable_analyses(self):
        path = Path(__file__).parent.parent / 'tests-nmr' / 'json-schema' / 'nmr-data-processing-report-schema-v5.json'
        schema = json.loads(path.read_text(encoding='utf-8'))
        Draft7Validator.check_schema(schema)
        stats = schema['definitions']['output_statistics']['properties']
        for key in ('entry_title', 'chem_shift_summary', 'restraint_summary'):
            validator = Draft7Validator(stats[key])
            with self.subTest(key=key):
                validator.validate(None)
                with self.assertRaises(ValidationError):
                    validator.validate([])
        title = Draft7Validator(stats['entry_title'])
        title.validate('Original entry title')
        with self.assertRaises(ValidationError):
            title.validate('')
        count = Draft7Validator(stats['assembly']['properties']['number_of_components'])
        count.validate(None)
        count.validate(1)
        for invalid in (0, -1, 'unknown'):
            with self.subTest(invalid=invalid), self.assertRaises(ValidationError):
                count.validate(invalid)
        name = Draft7Validator(stats['assembly']['properties']['name'])
        name.validate(None)
        name.validate('Original assembly')
        with self.assertRaises(ValidationError):
            name.validate(42)
