import React, { PureComponent } from "react";
import { connect } from "react-redux";
import { Drawer, H3 } from "@blueprintjs/core";

@connect((state) => {
  return { singleValueCategories: state.metadata.singleValueCategories };
})
class InfoDrawer extends PureComponent {
  constructor(props) {
    super(props);

    this.state = {
      isOpen: false,
    };
  }

  handleClick = () => {
    this.setState((state) => {
      return { isOpen: !state.isOpen };
    });
  };

  handleKeyPress = (e) => {
    if (e.key === "Enter") {
      this.handleCategoryClick();
    }
  };

  render() {
    const { isOpen } = this.state;
    const {
      position,
      aboutURL,
      datasetTitle,
      title,
      children,
      singleValueCategories,
    } = this.props;
    return (
      <div
        role="menuitem"
        tabIndex="0"
        onKeyPress={this.handleKeyPress}
        onClick={this.handleClick}
      >
        {children}
        <Drawer {...{ isOpen, position, title }}>
          <a
            style={{ width: 185 }}
            href={aboutURL}
            data-testid="header"
            target="_blank"
            rel="noopener noreferrer"
          >
            {datasetTitle}
          </a>
          <H3>Dataset Metadata</H3>
          {Array.from(singleValueCategories).map((pair) => {
            return (
              <div key={pair[0]} style={{ marginBottom: 10, marginTop: 4 }}>
                <span style={{ maxWidth: 150, fontWeight: 700 }}>
                  {pair[0]}
                </span>
                <span style={{ maxWidth: 150 }}>{`: ${pair[1]}`}</span>
              </div>
            );
          })}
        </Drawer>
      </div>
    );
  }
}
export default InfoDrawer;
