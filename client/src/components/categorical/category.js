import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import { FaChevronRight, FaChevronDown } from "react-icons/fa";
import { Button, Tooltip } from "@blueprintjs/core";

import * as globals from "../../globals";
import Value from "./value";
import alphabeticallySortedValues from "./util";

@connect(state => ({
  colorAccessor: state.controls.colorAccessor,
  categoricalSelectionState: state.controls.categoricalSelectionState
}))
class Category extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      isChecked: true,
      isExpanded: false
    };
  }

  componentDidUpdate() {
    const { categoricalSelectionState, metadataField } = this.props;
    const cat = categoricalSelectionState[metadataField];
    const categoryCount = {
      // total number of options in this category
      totalOptionCount: cat.numOptions,
      // number of selected options in this category
      selectedOptionCount: _.reduce(
        cat.optionSelected,
        (res, cond) => (cond ? res + 1 : res),
        0
      )
    };
    if (categoryCount.selectedOptionCount === categoryCount.totalOptionCount) {
      /* everything is on, so not indeterminate */
      this.checkbox.indeterminate = false;
    } else if (categoryCount.selectedOptionCount === 0) {
      /* nothing is on, so no */
      this.checkbox.indeterminate = false;
    } else if (
      categoryCount.selectedOptionCount < categoryCount.totalOptionCount
    ) {
      /* to be explicit... */
      this.checkbox.indeterminate = true;
    }
  }

  handleColorChange() {
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "color by categorical metadata",
      colorAccessor: metadataField
    });
  }

  toggleAll() {
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "categorical metadata filter all of these",
      metadataField
    });
    this.setState({ isChecked: true });
  }

  toggleNone() {
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "categorical metadata filter none of these",
      metadataField
    });
    this.setState({ isChecked: false });
  }

  handleToggleAllClick() {
    const { isChecked } = this.state;
    // || this.checkbox.indeterminate === false
    if (isChecked) {
      this.toggleNone();
    } else if (!isChecked) {
      this.toggleAll();
    }
  }

  renderCategoryItems() {
    const { categoricalSelectionState, metadataField } = this.props;

    const cat = categoricalSelectionState[metadataField];
    const optTuples = alphabeticallySortedValues([...cat.optionIndex]);
    return _.map(optTuples, (tuple, i) => (
      <Value
        key={tuple[1]}
        metadataField={metadataField}
        optionIndex={tuple[1]}
        i={i}
      />
    ));
  }

  render() {
    const { isExpanded, isChecked } = this.state;
    const {
      metadataField,
      colorAccessor,
      categoricalSelectionState
    } = this.props;
    const { isTruncated } = categoricalSelectionState[metadataField];
    return (
      <div
        style={{
          maxWidth: globals.maxControlsWidth
        }}
      >
        <div
          style={{
            display: "flex",
            justifyContent: "space-between",
            alignItems: "baseline"
          }}
        >
          <div
            style={{
              display: "flex",
              justifyContent: "flex-start",
              alignItems: "baseline"
            }}
          >
            <label className="bp3-control bp3-checkbox">
              <input
                onChange={this.handleToggleAllClick.bind(this)}
                ref={el => {
                  this.checkbox = el;
                  return el;
                }}
                checked={isChecked}
                type="checkbox"
              />
              <span className="bp3-control-indicator" />
              {""}
            </label>

            <span
              style={{
                cursor: "pointer",
                display: "inline-block"
              }}
              onClick={() => {
                this.setState({ isExpanded: !isExpanded });
              }}
            >
              {metadataField}
              {isExpanded ? (
                <FaChevronDown style={{ fontSize: 10, marginLeft: 5 }} />
              ) : (
                <FaChevronRight style={{ fontSize: 10, marginLeft: 5 }} />
              )}
            </span>
          </div>
          <Tooltip content="Use as color scale" position="bottom">
            <Button
              onClick={this.handleColorChange.bind(this)}
              active={colorAccessor === metadataField}
              intent={colorAccessor === metadataField ? "primary" : "none"}
              icon={"tint"}
            />
          </Tooltip>
        </div>
        <div style={{ marginLeft: 26 }}>
          {isExpanded ? this.renderCategoryItems() : null}
        </div>
        <div>
          {isExpanded && isTruncated ? (
            <p style={{ paddingLeft: 15 }}>... truncated list ...</p>
          ) : null}
        </div>
      </div>
    );
  }
}

export default Category;
