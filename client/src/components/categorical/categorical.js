// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { connect } from "react-redux";

import { FaChevronRight } from "react-icons/fa";
import { FaChevronDown } from "react-icons/fa";
import { FaPaintBrush } from "react-icons/fa";

import * as globals from "../../globals";
import Value from "./value";
import { alphabeticallySortedValues } from "./util";

@connect(state => ({
  colorAccessor: state.controls.colorAccessor,
  categoricalAsBooleansMap: state.controls.categoricalAsBooleansMap
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
    const { categoricalAsBooleansMap, metadataField } = this.props;

    const valuesAsBool = _.values(categoricalAsBooleansMap[metadataField]);
    /* count categories toggled on by counting true values */
    const categoriesToggledOn = _.values(valuesAsBool).filter(v => v).length;

    if (categoriesToggledOn === valuesAsBool.length) {
      /* everything is on, so not indeterminate */
      this.checkbox.indeterminate = false;
    } else if (categoriesToggledOn === 0) {
      /* nothing is on, so no */
      this.checkbox.indeterminate = false;
    } else if (categoriesToggledOn < valuesAsBool.length) {
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
    const { dispatch, metadataField, value } = this.props;
    dispatch({
      type: "categorical metadata filter none of these",
      metadataField,
      value
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
    const { values, metadataField } = this.props;
    return _.map(alphabeticallySortedValues(values), (v, i) => (
      <Value
        key={v}
        metadataField={metadataField}
        count={values[v]}
        value={v}
        i={i}
      />
    ));
  }

  render() {
    const { isExpanded, isChecked } = this.state;
    const { metadataField, colorAccessor } = this.props;
    return (
      <div
        style={{
          // display: "flex",
          // alignItems: "baseline",
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
          <p
            style={{
              // flexShrink: 0,
              fontWeight: 500,
              // textAlign: "right",
              // fontFamily: globals.accentFont,
              // fontStyle: "italic",
              margin: "3px 10px 3px 0px"
            }}
          >
            <span
              style={{
                cursor: "pointer",
                display: "inline-block",
                position: "relative",
                top: -1
              }}
              onClick={() => {
                this.setState({ isExpanded: !isExpanded });
              }}
            >
              {isExpanded ? <FaChevronDown /> : <FaChevronRight />}
            </span>
            {metadataField}
            <input
              onChange={this.handleToggleAllClick.bind(this)}
              ref={el => (this.checkbox = el)}
              checked={isChecked}
              type="checkbox"
            />
            <span
              onClick={this.handleColorChange.bind(this)}
              style={{
                fontSize: 14,
                marginLeft: 4,
                // padding: this.props.colorAccessor === this.props.metadataField ? 3 : "auto",
                borderRadius: 3,
                color:
                  colorAccessor === metadataField
                    ? globals.brightBlue
                    : "black",
                // backgroundColor: this.props.colorAccessor === this.props.metadataField ? globals.brightBlue : "inherit",
                display: "inline-block",
                position: "relative",
                top: -1,
                cursor: "pointer"
              }}
            >
              <FaPaintBrush />
            </span>
          </p>
        </div>
        <div>{isExpanded ? this.renderCategoryItems() : null}</div>
      </div>
    );
  }
}

@connect(state => {
  const ranges = _.get(state.controls.world, "summary.obs", null);
  return {
    ranges
  };
})
class Categories extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  render() {
    const { ranges } = this.props;
    if (!ranges) return null;

    return (
      <div
        style={{
          width: 310,
          marginRight: 40,
          paddingRight: 20,
          flexShrink: 0
          // height: 700,
          // overflow: "auto",
        }}
      >
        {_.map(ranges, (value, key) => {
          const isColorField = key.includes("color") || key.includes("Color");
          if (value.options && !isColorField && key !== "name") {
            return (
              <Category key={key} metadataField={key} values={value.options} />
            );
          }
        })}
      </div>
    );
  }
}

export default Categories;
