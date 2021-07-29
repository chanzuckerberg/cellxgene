import React from "react";
import { InputGroup, MenuItem, Keys } from "@blueprintjs/core";
import { Suggest } from "@blueprintjs/select";
import fuzzysort from "fuzzysort";

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
type State = any;

// eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
export default class LabelInput extends React.PureComponent<{}, State> {
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

  // eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
  constructor(props: {}) {
    super(props);

    // @ts-expect-error ts-migrate(2339) FIXME: Property 'label' does not exist on type '{}'.
    const { label } = props;
    const query = label || "";
    const queryResults = this.filterLabels(query);
    this.state = {
      query,
      queryResults,
    };
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleQueryChange = (query: any, event: any) => {
    // https://github.com/palantir/blueprint/issues/2983
    if (!event) return;

    const queryResults = this.filterLabels(query);
    this.setState({
      query,
      queryResults,
    });

    // @ts-expect-error ts-migrate(2339) FIXME: Property 'onChange' does not exist on type 'Readon... Remove this comment to see the full error message
    const { onChange } = this.props;
    if (onChange) onChange(query, event);
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleItemSelect = (item: any, event: any) => {
    /* only report the select if not already reported via onChange() */
    const { target } = item;
    const { query } = this.state;
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'onSelect' does not exist on type 'Readon... Remove this comment to see the full error message
    const { onSelect } = this.props;
    if (target !== query && onSelect) onSelect(target, event);
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleKeyDown = (e: any) => {
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleChange = (e: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'onChange' does not exist on type 'Readon... Remove this comment to see the full error message
    const { onChange } = this.props;
    if (onChange) onChange(e.target.value);
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  renderLabelSuggestion = (
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    queryResult: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    { handleClick, modifiers }: any
  ) => {
    if (queryResult.newLabel) {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'newLabelMessage' does not exist on type ... Remove this comment to see the full error message
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  filterLabels(query: any) {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'labelSuggestions' does not exist on type... Remove this comment to see the full error message
    const { labelSuggestions } = this.props;
    if (!labelSuggestions) return [];

    /* empty query is wildcard */
    if (query === "") {
      return (
        labelSuggestions
          .slice(0, LabelInput.QueryResultLimit)
          // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          .map((l: any) => ({
            target: l,
            score: -10000,
          }))
      );
    }

    /* else, do a fuzzy query */
    const options = {
      limit: LabelInput.QueryResultLimit,
      threshold: -10000, // don't return bad results
    };
    let queryResults = fuzzysort.go(query, labelSuggestions, options);
    /* exact match will always be first in list */
    if (query !== "" && queryResults[0]?.target !== query)
      // @ts-expect-error ts-migrate(2322) FIXME: Type '{ target: any; newLabel: true; }' is not ass... Remove this comment to see the full error message
      queryResults = [{ target: query, newLabel: true }, ...queryResults];

    return queryResults;
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    const { props } = this;
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'labelSuggestions' does not exist on type... Remove this comment to see the full error message
    const { labelSuggestions, label, autoFocus = true } = props;
    const suggestEnabled = !!labelSuggestions && labelSuggestions.length > 0;

    if (!suggestEnabled) {
      return (
        <InputGroup
          autoFocus={autoFocus}
          // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          {...(props as any).inputProps} // eslint-disable-line react/jsx-props-no-spreading --- Allows for modularity
          value={label}
          onChange={this.handleChange}
        />
      );
    }

    const popoverProps = {
      minimal: true,
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      ...(props as any).popoverProps,
    };
    const inputProps = {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      ...(props as any).inputProps,
      autoFocus: false,
    };
    const { queryResults } = this.state;
    return (
      <>
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
          // @ts-expect-error ts-migrate(2322) FIXME: Type '{ fill: true; inputValueRenderer: (i: any) =... Remove this comment to see the full error message
          onKeyDown={this.handleKeyDown}
        />
      </>
    );
  }
}
