import React from "react";
import { InputGroup, MenuItem, Keys } from "@blueprintjs/core";
import { Suggest } from "@blueprintjs/select";
import fuzzysort from "fuzzysort";

export default class LabelInput extends React.PureComponent {
  /*
  Input widget for text labels, which acts like an InputGroup, but will also 
  accept a suggestion list (of labels), with sublime-like suggest search.

  Properties:

  * labelSuggestions -- array of suggested lables. May be array of string or
    other objects. If array of objects, specify `labelKey`. If null, suggestion
    mode is disabled.
  * onSelect -- optional, callback upon selection of item from labelSuggestions.
    (label) => void
  * onChange -- optional, callback upon change in text input.  (label) => void
  * label -- component value, for controlled use
  * newLabelMessage -- text to display when user enters a label not in labelSuggestions
    (not used if suggestion mode disabled)
  * inputProps -- will be passed to InputGroup
  * popoverProps -- will be passed to <Suggest>
  */

  /* maxinum number of suggestions */
  static QueryResultLimit = 100;

  constructor(props) {
    super(props);

    const { label } = props;
    const query = label || "";
    const queryResults = this.filterLabels(query);
    this.state = {
      query,
      queryResults,
    };
  }

  handleQueryChange = (query, event) => {
    // https://github.com/palantir/blueprint/issues/2983
    if (!event) return;

    const queryResults = this.filterLabels(query);
    this.setState({
      query,
      queryResults,
    });

    const { onChange } = this.props;
    if (onChange) onChange(query, event);
  };

  handleItemSelect = (item, event) => {
    /* only report the select if not already reported via onChange() */
    const { target } = item;
    const { query } = this.state;
    const { onSelect } = this.props;
    if (target !== query && onSelect) onSelect(target, event);
  };

  handleKeyDown = (e) => {
    /* 
    prevent these events from propagating to containing form/dialog
    and causing further side effects (eg, closing dialog, submitting
    form, etc).
    */
    const { keyCode } = e;
    if (keyCode === Keys.ENTER || keyCode === Keys.ESCAPE) {
      e.preventDefault();
    }
    if (keyCode === Keys.ESCAPE) {
      e.stopPropagation();
    }
  };

  handleChange = (e) => {
    const { onChange } = this.props;
    if (onChange) onChange(e.target.value);
  };

  renderLabelSuggestion = (queryResult, { handleClick, modifiers }) => {
    if (queryResult.newLabel) {
      const { newLabelMessage } = this.props;
      return (
        <MenuItem
          icon="flag"
          active={modifiers.active}
          disabled={modifiers.disabled}
          key={queryResult.target}
          onClick={handleClick}
          text={<em>{queryResult.target}</em>}
          label={newLabelMessage || "New label"}
        />
      );
    }
    return (
      <MenuItem
        active={modifiers.active}
        disabled={modifiers.disabled}
        key={queryResult.target}
        onClick={handleClick}
        text={queryResult.target}
      />
    );
  };

  filterLabels(query) {
    const { labelSuggestions } = this.props;
    if (!labelSuggestions) return [];

    /* empty query is wildcard */
    if (query === "") {
      return labelSuggestions
        .slice(0, LabelInput.QueryResultLimit)
        .map((l) => ({
          target: l,
          score: -10000,
        }));
    }

    /* else, do a fuzzy query */
    const options = {
      limit: LabelInput.QueryResultLimit,
      threshold: -10000, // don't return bad results
    };
    let queryResults = fuzzysort.go(query, labelSuggestions, options);
    /* exact match will always be first in list */
    if (query !== "" && queryResults[0]?.target !== query)
      queryResults = [{ target: query, newLabel: true }, ...queryResults];

    return queryResults;
  }

  render() {
    const { props } = this;
    const { labelSuggestions, label, autoFocus = true } = props;
    const suggestEnabled = !!labelSuggestions && labelSuggestions.length > 0;

    if (!suggestEnabled) {
      return (
        <InputGroup
          autoFocus={autoFocus}
          {...props.inputProps} // eslint-disable-line react/jsx-props-no-spreading --- Allows for modularity
          value={label}
          onChange={this.handleChange}
        />
      );
    }

    const popoverProps = {
      minimal: true,
      ...props.popoverProps,
    };
    const inputProps = {
      ...props.inputProps,
      autoFocus: false,
    };
    const { queryResults } = this.state;
    return (
      <Suggest
          fill
          inputValueRenderer={(i) => i.target}
          items={queryResults}
          itemRenderer={this.renderLabelSuggestion}
          onItemSelect={this.handleItemSelect}
          query={label}
          onQueryChange={this.handleQueryChange}
          popoverProps={popoverProps}
          inputProps={inputProps}
          onKeyDown={this.handleKeyDown}
        />
    );
  }
}
