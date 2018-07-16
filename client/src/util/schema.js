// jshint esversion: 6

// In the case where the REST server does not implement data schema
// declaration, we attempt to deduce it by sniffing the data.
//
export function createSchemaByDataSniffing(ranges) {
  let schema = {};
  _.forEach(ranges, (value, key) => {
    schema[key] = {
      displayname: key,
      variabletype: value.options ? "categorical" : "continuous"
    };

    // Metadata field type is inferred by sniffing the data. This has some risks.
    // Caveats:
    //  * Values have been converted to native JS objects by the JSON parser.
    //  * Lots of assumptions about he REST API behaving properly (eg, min/max
    //    are the same type, etc).
    let type;
    if (schema[key].variabletype === "continuous" && value.range) {
      // Use min/max as a proxy for all data.
      const min = value.range.min;
      const max = value.range.max;
      type =
        typeof min !== "number" || typeof max !== "number"
          ? "string"
          : Number.isSafeInteger(min) && Number.isSafeInteger(max)
            ? "int"
            : "float";
    } else {
      // use an option value as a proxy for all data
      const aVal = value.options[0];
      type =
        typeof aVal !== "number"
          ? "string"
          : Number.isSafeInteger(aVal) ? "int" : "float";
    }
    schema[key].type = type;
  });
  return schema;
}
